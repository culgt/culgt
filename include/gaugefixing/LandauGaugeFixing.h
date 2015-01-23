/**
 * LandauGaugeFixing.h
 *
 *  Created on: Mar 19, 2014
 *      Author: vogt
 */

#ifndef LANDAUGAUGEFIXING_H_
#define LANDAUGAUGEFIXING_H_

#include <assert.h>
#include "../lattice/parameterization_types/SU2Vector4.h"
#include "../lattice/SubgroupIterator.h"
#include "../common/culgt_typedefs.h"
#include "../cuLGT1legacy/Chronotimer.h"
#include <iostream>
#include "../cuLGT1legacy/Reduction.hxx"
#include "../cudacommon/cuda_error.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "algorithms/OrUpdate.h"
#include "algorithms/SaUpdate.h"
#include "algorithms/MicroUpdate.h"
#include "gaugetypes/LandauCoulombGaugeType.h"
#include "RunInfo.h"
#include "GaugeStats.h"
#include <string>

#include "gaugefixing_thread_types.h"
#include "../lattice/GaugeConfigurationCudaHelper.h"
#include "../util/performance/TunableObject.h"
#include "FullGaugeFixingOverrelaxation.h"
#include "FullGaugeFixingSimulatedAnnealing.h"
#include "GaugeFixingSaOr.h"
#include "RandomGaugeTrafo.h"
#include "../lattice/GlobalLink.h"

using std::string;

namespace culgt
{

namespace LandauGaugefixingKernel
{

template<typename GlobalLinkType, typename LocalLinkType>  __global__ void generateGaugeQualityPerSite( typename GlobalLinkType::PARAMTYPE::TYPE* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, lat_index_t* nn, double *dGff, double *dA )
{
	typename GlobalLinkType::PATTERNTYPE::SITETYPE site( dim, nn );

	lat_index_t index = blockIdx.x * blockDim.x + threadIdx.x;

	LocalLinkType Sum;
	Sum.zero();

	double gff = 0;
	for( int mu = 0; mu < 4; mu++ )
	{
		LocalLinkType temp;

		site.setLatticeIndex( index );
		GlobalLinkType globalLink( U, site, mu );
		temp = globalLink;
		Sum += temp;
		gff += temp.reTrace();

		site.setNeighbour(mu,false);
		GlobalLinkType globDw( U, site, mu );
		temp = globDw;
		Sum -= temp;
	}

	Sum -= Sum.trace()/(LocalLinkType::PARAMTYPE::REALTYPE)(GlobalLinkType::PARAMTYPE::NC);

	LocalLinkType SumHerm;
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	dA[index] = Sum.normFrobeniusSquared();
	dGff[index] = gff;
}

}

template<typename PatternType, typename LocalLinkType> class LandauGaugefixing: public GaugeFixingSaOr
{
public:
	typedef GlobalLink<PatternType,true> GlobalLinkType;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	LandauGaugefixing( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, long seed ) : GaugeFixingSaOr( dim.getSize() ), dim(dim), U(U), overrelaxation( &this->U, dim, 1.5, seed ), simulatedAnnealing( &this->U, dim, seed, 1.0 ), microcanonical( &this->U, dim, seed ), seed(seed)
	{
		// TODO we should assert here that we use GPUPatternParitySplit!
	}

	~LandauGaugefixing()
	{
	}

	GaugeStats getGaugeStats( GaugeFieldDefinition defintion = GAUGEFIELD_STANDARD )
	{
		GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupNoSplit( dim, false );
		LandauGaugefixingKernel::generateGaugeQualityPerSite<GlobalLinkType,LocalLinkType><<<setupNoSplit.getGridSize(),setupNoSplit.getBlockSize()>>>( U, dim, SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), dGff, dA );
		CUDA_LAST_ERROR( "generateGaugeQualityPerSite" );

		Reduction<double> reducer(dim.getSize());
		double dAAvg = reducer.reduceAll( dA )/(double)dim.getSize()/(double)(GlobalLinkType::PATTERNTYPE::PARAMTYPE::NC);
		double dGffAvg = reducer.reduceAll( dGff )/(double)dim.getSize()/(double)(GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim*GlobalLinkType::PATTERNTYPE::PARAMTYPE::NC);

		return GaugeStats( dGffAvg, dAAvg );
	}


	void runOverrelaxation( float orParameter  )
	{
		overrelaxation.setOrParameter( orParameter );
		overrelaxation.run();
	}

	RunInfo getRunInfoOverrelaxation( double time, long iter )
	{
		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<LANDAU> >( dim.getSize(), time, iter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	void runMicrocanonical()
	{
		microcanonical.run();
	}

	void runSimulatedAnnealing( float temperature )
	{
		simulatedAnnealing.setTemperature( temperature );
		simulatedAnnealing.run();
	}

	RunInfo getRunInfoSimulatedAnnealing( double time, long iterSa, long iterMicro )
	{
		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<LANDAU> >( dim.getSize(), time, iterSa, SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops, iterMicro, MicroUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	template<typename RNG> void orstepsAutoTune( float orParameter = 1.5, int iter = 1000 )
	{
		overrelaxation.setOrParameter( orParameter );
		overrelaxation.tune( iter );
	}

	template<typename RNG> void sastepsAutoTune( float temperature = 1.0, int iter = 1000 )
	{
		simulatedAnnealing.setTemperature( temperature );
		simulatedAnnealing.tune( iter );
	}

	template<typename RNG> void microcanonicalAutoTune( int iter = 1000 )
	{
		microcanonical.tune( iter );
	}

	void randomTrafo()
	{
		RandomGaugeTrafo<PatternType,LocalLinkType>::randomTrafo( U, dim, seed );
	}

	void reproject()
	{
		GaugeConfigurationCudaHelper<T>::template reproject<typename GlobalLinkType::PATTERNTYPE,LocalLinkType, PhiloxWrapper<REALT> >( U, dim, seed, PhiloxWrapper<REALT>::getNextCounter() );
	}

	void allocateCopyMemory()
	{
		GaugeConfigurationCudaHelper<T>::allocateMemory( &UBest, GlobalLinkType::getArraySize(dim) );
		GaugeConfigurationCudaHelper<T>::allocateMemory( &UClean, GlobalLinkType::getArraySize(dim) );
	}

	void freeCopyMemory()
	{
		GaugeConfigurationCudaHelper<T>::freeMemory( UBest );
		GaugeConfigurationCudaHelper<T>::freeMemory( UClean );
	}

	void saveCopy()
	{
		CUDA_SAFE_CALL( cudaMemcpy( UBest, U, GlobalLinkType::getArraySize(dim)*sizeof(T), cudaMemcpyDeviceToDevice ), "cudaMemcpy in saveCopy (device to device)" );
	}

	void writeBackCopy()
	{
		CUDA_SAFE_CALL( cudaMemcpy( U, UBest, GlobalLinkType::getArraySize(dim)*sizeof(T), cudaMemcpyDeviceToDevice ), "cudaMemcpy in saveCopy (device to device)" );
	}

	void storeCleanCopy()
	{
		CUDA_SAFE_CALL( cudaMemcpy( UClean, U, GlobalLinkType::getArraySize(dim)*sizeof(T), cudaMemcpyDeviceToDevice ), "cudaMemcpy in saveCopy (device to device)" );
	}

	void takeCleanCopy()
	{
		CUDA_SAFE_CALL( cudaMemcpy( U, UClean, GlobalLinkType::getArraySize(dim)*sizeof(T), cudaMemcpyDeviceToDevice ), "cudaMemcpy in saveCopy (device to device)" );
	}

private:
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim;

	T* U;

	T* UBest;
	T* UClean;

	FullGaugeFixingOverrelaxation<PatternType,LocalLinkType,LandauCoulombGaugeType<LANDAU>,false > overrelaxation;
	FullGaugeFixingOverrelaxation<PatternType,LocalLinkType,LandauCoulombGaugeType<LANDAU>,true > microcanonical;
	FullGaugeFixingSimulatedAnnealing<PatternType,LocalLinkType,LandauCoulombGaugeType<LANDAU> > simulatedAnnealing;

	long seed;
};


}

#endif /* LANDAUGAUGEFIXING_H_ */
