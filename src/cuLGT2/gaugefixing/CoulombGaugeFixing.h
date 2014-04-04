/**
 *  Created on: Mar 19, 2014
 *      Author: vogt
 */

#ifndef COULOMBGAUGEFIXING_H_
#define COULOMBGAUGEFIXING_H_
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

using std::string;

namespace culgt
{

//enum GaugeFixingThreadsPerSite{ SINGLE_THREAD_PER_SITE, FOUR_THREAD_PER_SITE, EIGHT_THREAD_PER_SITE, TIMESLICE };


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
	Sum -= Sum.trace()/Real(GlobalLinkType::PARAMTYPE::NC);

	LocalLinkType SumHerm;
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	dA[index] = Sum.normFrobeniusSquared();
	dGff[index] = gff;
}


template<typename GlobalLinkType, typename GlobalLinkType2, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(8*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep8( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter)
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing8Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
	gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename GlobalLinkType2,  typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(4*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing4Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
	gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename GlobalLinkType2, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(8*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelSaStep8( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE temperature, int rngSeed, int rngCounter)
{
	PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	typedef SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm saupdate( temperature, rng );
	GaugeFixing8Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
	gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename GlobalLinkType2,  typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(4*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelSaStep4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE temperature, int rngSeed, int rngCounter )
{
	PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	typedef SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm saupdate( temperature, rng );
	GaugeFixing4Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
	gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
}

}

template<typename GlobalLinkType, typename LocalLinkType> class CoulombGaugefixing
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );

	CoulombGaugefixing( T* Ut, T* UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim ) : dimTimeslice(dim), totalTime(0), totalIter(0)
	{
		CUDA_SAFE_CALL( cudaMalloc( &dA, dim.getSize()*sizeof(double) ), "malloc dA");
		CUDA_SAFE_CALL( cudaMalloc( &dGff, dim.getSize()*sizeof(double) ), "malloc dGff");

		setTimeslice( Ut, UtDown );
	}

	~CoulombGaugefixing()
	{
		cudaFree( dA );
		cudaFree( dGff );
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

	template<GaugeFixingThreadsPerSite ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> RunInfo orsteps( float orParameter, int iter )
	{
		if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
		{
			cudaFuncSetCacheConfig( CoulombGaugefixingKernel::kernelOrStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
		{
			cudaFuncSetCacheConfig( CoulombGaugefixingKernel::kernelOrStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dimTimeslice, true, SitesPerBlock );

		Chronotimer timer;
		timer.reset();
		timer.start();
		cudaDeviceSynchronize();

		for( int i = 0; i < iter; i++ )
		{
			if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
			{
				CoulombGaugefixingKernel::kernelOrStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, orParameter );
				CoulombGaugefixingKernel::kernelOrStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, orParameter );
			}
			else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
			{
				CoulombGaugefixingKernel::kernelOrStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, orParameter );
				CoulombGaugefixingKernel::kernelOrStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, orParameter );
			}
			else
			{
				assert( false );
			}
		}

		cudaDeviceSynchronize();
		timer.stop();
		CUDA_LAST_ERROR( "kernelOrStep" );

		totalTime += timer.getTime();
		totalIter += iter;

		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<COULOMB> >( dimTimeslice.getSize(), timer.getTime(), iter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	template<GaugeFixingThreadsPerSite ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> RunInfo sasteps( float saMax, float saMin, int steps, int seed )
	{
		if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
		{
			cudaFuncSetCacheConfig( CoulombGaugefixingKernel::kernelSaStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
		{
			cudaFuncSetCacheConfig( CoulombGaugefixingKernel::kernelSaStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dimTimeslice, true, SitesPerBlock );

		Chronotimer timer;
		timer.reset();
		timer.start();
		cudaDeviceSynchronize();

		float tStep = (saMax-saMin)/(float)steps;
		float temperature = saMax;
		for( int i = 0; i < steps; i++ )
		{
			if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
			{
				CoulombGaugefixingKernel::kernelSaStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, temperature, seed, PhiloxWrapper<REALT>::getNextCounter() );
				CoulombGaugefixingKernel::kernelSaStep8<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, temperature, seed, PhiloxWrapper<REALT>::getNextCounter() );
			}
			else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
			{
				CoulombGaugefixingKernel::kernelSaStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, temperature, seed, PhiloxWrapper<REALT>::getNextCounter() );
				CoulombGaugefixingKernel::kernelSaStep4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, temperature, seed, PhiloxWrapper<REALT>::getNextCounter() );
			}
			else
			{
				assert( false );
			}
			temperature -= tStep;
		}

		cudaDeviceSynchronize();
		timer.stop();
		CUDA_LAST_ERROR( "kernelOrStep" );

		totalTime += timer.getTime();
		totalIter += steps;

		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<COULOMB> >( dimTimeslice.getSize(), timer.getTime(), steps, SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

private:
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice;
	T* Ut;
	T* UtDown;
	double* dA;
	double* dGff;

	double totalTime;
	int totalIter;
};


}

#endif /* COULOMBGAUGEFIXING_H_ */
