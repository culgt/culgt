/**
 *  Created on: Mar 19, 2014
 *      Author: vogt
 */

#ifndef COULOMBGAUGEFIXING_H_
#define COULOMBGAUGEFIXING_H_
#include "gaugefixing_thread_types.h"
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
#include "gaugetypes/LandauCoulombGaugeType.h"
#include "RunInfo.h"
#include "GaugeStats.h"
#include <string>
#include "../lattice/GlobalLink.h"
#include "../util/rng/PhiloxWrapper.h"

#include "GaugeSettings.h"

#include "CoulombGaugeFixingOverrelaxation.h"
#include "CoulombGaugeFixingSimulatedAnnealing.h"
#include "AutoTuner.h"
#include "RandomGaugeTrafo.h"

#include "GaugeFixingSaOr.h"

using std::string;


namespace culgt
{



namespace CoulombGaugefixingKernel
{

template<typename GlobalLinkType, typename LocalLinkType>  __global__ void generateGaugeQualityPerSite( typename GlobalLinkType::PARAMTYPE::TYPE* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, lat_index_t* nn, double *dGff, double *dA )
{
	typename GlobalLinkType::PATTERNTYPE::SITETYPE site( dim, nn );

	lat_index_t index = blockIdx.x * blockDim.x + threadIdx.x;

	LocalLinkType Sum;
	Sum.zero();

	double gff = 0;
	for( int mu = 1; mu < 4; mu++ )
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

	// TODO: verify that the following statement indeed drops out, then remove it
	Sum -= Sum.trace()/(LocalLinkType::PARAMTYPE::REALTYPE)(GlobalLinkType::PARAMTYPE::NC);

	LocalLinkType SumHerm;
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	dA[index] = Sum.normFrobeniusSquared();
	dGff[index] = gff;
}

}

template<typename GlobalLinkType, typename LocalLinkType> class CoulombGaugefixing: public GaugeFixingSaOr
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );

	CoulombGaugefixing( T* Ut, T* UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, long seed ) : GaugeFixingSaOr( dim.getSize() ), dimTimeslice(dim), totalTime(0), totalIter(0), overrelaxation( Ut, UtDown, dim, seed, 1.5 ), simulatedAnnealing( Ut, UtDown, dim, seed, 1. ), seed(seed)
	{
		setTimeslice( Ut, UtDown );
	}

	void setTimeslice( T* Ut, T* UtDown )
	{
		this->Ut = Ut;
		GlobalLinkType::bindTexture( Ut, GlobalLinkType::getArraySize( dimTimeslice ) );
		this->UtDown = UtDown;
		GlobalLinkType2::bindTexture( UtDown, GlobalLinkType2::getArraySize( dimTimeslice ) );
	}

	GaugeStats getGaugeStats()
	{
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupNoSplit( dimTimeslice, false );
		CoulombGaugefixingKernel::generateGaugeQualityPerSite<GlobalLinkType,LocalLinkType><<<setupNoSplit.getGridSize(),setupNoSplit.getBlockSize()>>>( Ut, dimTimeslice, SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), dGff, dA );
		CUDA_LAST_ERROR( "generateGaugeQualityPerSite" );

		Reduction<double> reducer(dimTimeslice.getSize());
		double dAAvg = reducer.reduceAll( dA )/(double)dimTimeslice.getSize()/(double)(GlobalLinkType::PATTERNTYPE::PARAMTYPE::NC);
		double dGffAvg = reducer.reduceAll( dGff )/(double)dimTimeslice.getSize()/(double)((GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim-1)*GlobalLinkType::PATTERNTYPE::PARAMTYPE::NC);

		return GaugeStats( dGffAvg, dAAvg );
	}

	RunInfo getTotalRunInfo()
	{
		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<COULOMB> >( dimTimeslice.getSize(), totalTime, totalIter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	RunInfo runOverrelaxation( float orParameter, int iter, int id = -1 )
	{
		overrelaxation.setOrParameter( orParameter );

		overrelaxation.startTime();

		for( int i = 0; i < iter; i++ )
		{
			overrelaxation.run( id );
		}
		overrelaxation.stopTime();

		CUDA_LAST_ERROR( "kernelOrStep" );

		totalTime += overrelaxation.getTime();
		totalIter += iter;

		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<COULOMB> >( dimTimeslice.getSize(), overrelaxation.getTime(), iter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	template<typename RNG> void orstepsAutoTune( float orParameter = 1.5, int iter = 1000 )
	{
		overrelaxation.setOrParameter( orParameter );
		overrelaxation.tune();
	}

	RunInfo runSimulatedAnnealing( float saMax, float saMin, int steps, int id = -1 )
	{
		simulatedAnnealing.startTime();

		float tStep = (saMax-saMin)/(float)steps;
		float temperature = saMax;
		for( int i = 0; i < steps; i++ )
		{
			simulatedAnnealing.setTemperature( temperature );
			simulatedAnnealing.run( id );
			temperature -= tStep;
		}

		simulatedAnnealing.stopTime();
		CUDA_LAST_ERROR( "kernelSaStep" );

		totalTime += simulatedAnnealing.getTime();
		totalIter += steps;

		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<COULOMB> >( dimTimeslice.getSize(), simulatedAnnealing.getTime(), steps, SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	template<typename RNG> void sastepsAutoTune( float temperature = 1.0, int Id = -1 )
	{
		simulatedAnnealing.setTemperature( temperature );
		simulatedAnnealing.tune();
	}

	void randomTrafo()
	{
		RandomGaugeTrafo<GlobalLinkType,LocalLinkType>::randomTrafo( Ut, UtDown, dimTimeslice, seed );
	}



private:
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice;
	T* Ut;
	T* UtDown;


	CoulombGaugeFixingOverrelaxation<GlobalLinkType,LocalLinkType> overrelaxation;
	CoulombGaugeFixingSimulatedAnnealing<GlobalLinkType,LocalLinkType> simulatedAnnealing;

	double totalTime;
	int totalIter;

	int orOptimalTunedId;
	int saOptimalTunedId;

	long seed;
};


}

#endif /* COULOMBGAUGEFIXING_H_ */
