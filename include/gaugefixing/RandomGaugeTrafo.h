/**
 *  Created on: Mar 19, 2014
 *      Author: vogt
 */

#ifndef RANDOMGAUGETRAFO_H_
#define RANDOMGAUGETRAFO_H_
#include "../lattice/GlobalLink.h"
#include "../util/rng/PhiloxWrapper.h"
#include "gaugefixing_thread_types.h"
#include "algorithms/RandomUpdate.h"
#include "gaugetypes/RandomTrafoType.h"
#include "GaugeFixing4Threads.h"
#include "../cudacommon/cuda_error.h"
#include "../lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"

namespace culgt
{



namespace RandomGaugeTrafoKernel
{

template<typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock> __global__  void kernelRandom4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	typedef RandomUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm randomupdate( rng );
	GaugeFixing4Threads<Algorithm, RandomTrafoType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( randomupdate );
	gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
}

template<typename GlobalLinkType, typename GlobalLinkType2, typename LocalLinkType, int SitesPerBlock> __global__  void kernelRandom4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	typedef RandomUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm randomupdate( rng );
	GaugeFixing4Threads<Algorithm, RandomTrafoType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( randomupdate );
	gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
}

}


template<typename PatternType, typename LocalLinkType> class RandomGaugeTrafo
{
};

template<typename SiteType, typename ParamType, typename LocalLinkType> class RandomGaugeTrafo<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType >
{
public:
	typedef GlobalLink<GPUPatternParityPriority<SiteType,ParamType> > GlobalLinkType;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	static void randomTrafo( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, int seed )
	{
		const int SitesPerBlock = 128;

		GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
		CUDA_LAST_ERROR( "kernelRandomTrafo" );
	}

	static void randomTrafo( T* Ut, T* UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, int seed )
	{
		const int SitesPerBlock = 128;
		COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );

		GlobalLinkType::bindTexture( Ut, GlobalLinkType::getArraySize( dim ) );
		GlobalLinkType2::bindTexture( UtDown, GlobalLinkType2::getArraySize( dim ) );

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,GlobalLinkType2,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
		CUDA_LAST_ERROR( "kernelRandomTrafo" );
	}
};

template<typename SiteType, typename ParamType, typename LocalLinkType> class RandomGaugeTrafo<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType >
{
public:
	typedef GlobalLink<GPUPatternTimesliceParityPriority<SiteType,ParamType> > GlobalLinkType;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	typedef SiteIndex<SiteType::Ndim, FULL_SPLIT> SiteTypeTimeslice;
	typedef GlobalLink<GPUPatternParityPriority<SiteTypeTimeslice,ParamType> > GlobalLinkTypeTimeslice;

	static void randomTrafo( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, int seed )
	{
		COPY_GLOBALLINKTYPE( GlobalLinkTypeTimeslice, GlobalLinkTypeTimeslice2, 1 );
		const int SitesPerBlock = 128;

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim.getDimensionTimeslice(), true, SitesPerBlock );

		int timesliceArraySize = GlobalLinkType::getArraySize( dim.getDimensionTimeslice() );
		for( int t = 0; t < dim.getDimension(0); t++ )
		{
			int tdown = (t==0)?( dim.getDimension(0)-1):(t-1);

			GlobalLinkTypeTimeslice::bindTexture( &U[t*timesliceArraySize], timesliceArraySize );
			GlobalLinkTypeTimeslice2::bindTexture( &U[tdown*timesliceArraySize], timesliceArraySize );

			RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*timesliceArraySize], &U[tdown*timesliceArraySize], dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( dim.getDimensionTimeslice() ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
			RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( &U[t*timesliceArraySize], &U[tdown*timesliceArraySize], dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( dim.getDimensionTimeslice() ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
			CUDA_LAST_ERROR( "kernelRandomTrafo" );
		}

	}
};


}

#endif /* COULOMBGAUGEFIXING_H_ */
