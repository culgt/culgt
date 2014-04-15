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

}

template<typename GlobalLinkType, typename LocalLinkType> class RandomGaugeTrafo
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	RandomGaugeTrafo( T* U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim ) : dim(dim), U(U)
	{
	}

	~RandomGaugeTrafo()
	{
	}

	void randomTrafo( int seed )
	{
		const int SitesPerBlock = 32;

		GlobalLinkType::bindTexture( U, GlobalLinkType::getArraySize( dim ) );
//		cudaFuncSetCacheConfig( RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock>, cudaFuncCachePreferL1 );

		std::cout << "random trafo" << std::endl;

		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( dim, true, SitesPerBlock );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), false, seed, PhiloxWrapper<REALT>::getNextCounter() );
		RandomGaugeTrafoKernel::kernelRandom4<GlobalLinkType,LocalLinkType,SitesPerBlock><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*4,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( U, dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( dim ), true, seed, PhiloxWrapper<REALT>::getNextCounter() );
		CUDA_LAST_ERROR( "kernelRandomTrafo" );
	}

private:
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim;
	T* U;
};


}

#endif /* COULOMBGAUGEFIXING_H_ */
