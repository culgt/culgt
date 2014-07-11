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

template<typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock> __global__  void kernelRandom4( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	typedef RandomUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
	Algorithm randomupdate( rng );
	GaugeFixing4Threads<Algorithm, RandomTrafoType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( randomupdate );
	gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
}

}

template<typename GlobalLinkType, typename LocalLinkType> class RandomGaugeTrafo
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	static void randomTrafo( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, int seed )
	{
		const int SitesPerBlock = 32;

		GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );
//		cudaFuncSetCacheConfig( RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock>, cudaFuncCachePreferL1 );

//		std::cout << "random trafo" << std::endl;

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
		CUDA_LAST_ERROR( "kernelRandomTrafo" );
	}

	static void randomTrafo( T* Ut, T* UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, int seed )
	{
		const int SitesPerBlock = 32;
		COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );

		GlobalLinkType::bindTexture( Ut, GlobalLinkType::getArraySize( dim ) );
		GlobalLinkType2::bindTexture( UtDown, GlobalLinkType2::getArraySize( dim ) );

//		cudaFuncSetCacheConfig( RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock>, cudaFuncCachePreferL1 );

//		std::cout << "random trafo" << std::endl;

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( Ut, UtDown, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
		CUDA_LAST_ERROR( "kernelRandomTrafo" );
	}
};


}

#endif /* COULOMBGAUGEFIXING_H_ */
