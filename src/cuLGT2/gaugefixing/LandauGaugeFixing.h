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
#include "gaugetypes/LandauCoulombGaugeType.h"
#include "RunInfo.h"
#include "GaugeStats.h"
#include <string>

#include "gaugefixing_thread_types.h"


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

	// TODO: verify that the following statement indeed drops out, then remove it
	Sum -= Sum.trace()/(LocalLinkType::PARAMTYPE::REALTYPE)(GlobalLinkType::PARAMTYPE::NC);

	LocalLinkType SumHerm;
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	dA[index] = Sum.normFrobeniusSquared();
	dGff[index] = gff;
}


template<typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(8*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep8( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter)
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing8Threads<Algorithm, LandauCoulombGaugeType<LANDAU>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
	gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(4*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing4Threads<Algorithm, LandauCoulombGaugeType<LANDAU>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
	gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(4*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep4Timeslice( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Uup, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Udown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing4Threads<Algorithm, LandauCoulombGaugeType<LANDAU>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
	gaugefixing.applyAlgorithmTimeslice( Uup, Udown, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename GlobalLinkType2, typename LocalLinkType, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(8*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep8Timeslice( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Uup, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Udown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
{
	typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm orupdate( orParameter );
	GaugeFixing8Threads<Algorithm, LandauCoulombGaugeType<LANDAU>, GlobalLinkType, LocalLinkType, SitesPerBlock,GlobalLinkType2> gaugefixing( orupdate );
	gaugefixing.applyAlgorithmTimeslice( Uup, Udown, nn, latticeSize, parity );
}

}

template<typename GlobalLinkType, typename LocalLinkType> class LandauGaugefixing
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	LandauGaugefixing( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim ) : dim(dim), U(U), totalTime(0), totalIter(0)
	{
		// TODO we should assert here that we use GPUPatternParitySplit!
		CUDA_SAFE_CALL( cudaMalloc( &dA, dim.getSize()*sizeof(double) ), "malloc dA");
		CUDA_SAFE_CALL( cudaMalloc( &dGff, dim.getSize()*sizeof(double) ), "malloc dGff");
	}

	~LandauGaugefixing()
	{
		cudaFree( dA );
		cudaFree( dGff );
	}

	GaugeStats getGaugeStats()
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

	RunInfo getTotalRunInfo()
	{
		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<LANDAU> >( dim.getSize(), totalTime, totalIter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

	template<GaugeFixingThreadsPerSite ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> RunInfo orsteps( float orParameter, int iter )
	{

		if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
		{
			GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );
			cudaFuncSetCacheConfig( LandauGaugefixingKernel::kernelOrStep8<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
		{
			GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );
			cudaFuncSetCacheConfig( LandauGaugefixingKernel::kernelOrStep4<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
		}
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );

		Chronotimer timer;
		timer.reset();
		timer.start();
		cudaDeviceSynchronize();

		for( int i = 0; i < iter; i++ )
		{
			if( ThreadsPerSite == EIGHT_THREAD_PER_SITE )
			{
				LandauGaugefixingKernel::kernelOrStep8<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, orParameter );
				LandauGaugefixingKernel::kernelOrStep8<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, orParameter );
			}
			else if( ThreadsPerSite == FOUR_THREAD_PER_SITE )
			{
				LandauGaugefixingKernel::kernelOrStep4<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, orParameter );
				LandauGaugefixingKernel::kernelOrStep4<GlobalLinkType,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, orParameter );
			}
//			else if( ThreadsPerSite == TIMESLICE )
//			{
//				/*
//				 * This is just a test part to prepare for Coulomb gauge (or maybe for multi-GPU)
//				 * When you try this use GPUPatternTimesliceParitySplit (with a Site that is TIMESLICE_SPLIT) as the layout of U (in the kernel the work is done on a timeslice U[t] that then has size (1,x,x,x) with GPUPatternParitySplit layout and the Site to be used in the kernel is FULL_SPLIT
//				 * Therefore, the neighbour table needs to be FULL_SPLIT on a (1,x,x,x) lattice (see below).
//				 */
//				typedef GlobalLink<typename GlobalLinkType::PATTERNTYPE::TIMESLICE_PATTERNTYPE, GlobalLinkType::USETEXTURE> GlobalLinkTypeTimeslice;
////				typedef GlobalLink<typename GlobalLinkType::PATTERNTYPE::TIMESLICE_PATTERNTYPE, GlobalLinkType::USETEXTURE, 1> GlobalLinkTypeTimeslice2;
//
//
//				COPY_GLOBALLINKTYPE( GlobalLinkTypeTimeslice, GlobalLinkTypeTimeslice2, 1 );
//
//				LatticeDimension<4> dimTimeslice( 1, dim.getDimension(1), dim.getDimension(2), dim.getDimension(3) );
//				KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupTimeslice( dimTimeslice, true, SitesPerBlock );
//
//				lat_array_index_t arraySizeTimeslice = GlobalLinkTypeTimeslice::getArraySize( dimTimeslice );
//
//				for( int t = 0; t < dim.getDimension(0); t++ )
//				{
//					int tDown = (t>0)?(t-1):(dim.getDimension(0)-1);
//
//					/*
//					 * Binding textures in every step is surely not a good choice. But it is not necessary in Coulomb gauge (only once per timeslice)
//					 */
//					GlobalLinkTypeTimeslice::bindTexture( &U[t*arraySizeTimeslice], arraySizeTimeslice );
//					GlobalLinkTypeTimeslice2::bindTexture( &U[tDown*arraySizeTimeslice], arraySizeTimeslice );
////					LandauGaugefixingKernel::kernelOrStep4Timeslice<GlobalLinkTypeTimeslice,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupTimeslice.getGridSize(),setupTimeslice.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*arraySizeTimeslice], &U[tDown*arraySizeTimeslice], dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeTimeslice::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, orParameter );
////					LandauGaugefixingKernel::kernelOrStep4Timeslice<GlobalLinkTypeTimeslice,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupTimeslice.getGridSize(),setupTimeslice.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*arraySizeTimeslice], &U[tDown*arraySizeTimeslice], dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeTimeslice::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, orParameter );
//					LandauGaugefixingKernel::kernelOrStep8Timeslice<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupTimeslice.getGridSize(),setupTimeslice.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*arraySizeTimeslice], &U[tDown*arraySizeTimeslice], dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeTimeslice::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), false, orParameter );
//					LandauGaugefixingKernel::kernelOrStep8Timeslice<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupTimeslice.getGridSize(),setupTimeslice.getBlockSize()*8,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*arraySizeTimeslice], &U[tDown*arraySizeTimeslice], dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeTimeslice::PATTERNTYPE::SITETYPE>::getDevicePointer( dimTimeslice ), true, orParameter );
//				}
//			}
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

		return RunInfo::makeRunInfo<GlobalLinkType,LocalLinkType,LandauCoulombGaugeType<LANDAU> >( dim.getSize(), timer.getTime(), iter, OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE>::Flops );
	}

private:
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim;
	T* U;
	double* dA;
	double* dGff;

	double totalTime;
	int totalIter;
};


}

#endif /* LANDAUGAUGEFIXING_H_ */
