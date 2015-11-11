/*
 */

#ifndef FULLGAUGEFIXINGSIMULATEDANNEALINGKERNEL_H_
#define FULLGAUGEFIXINGSIMULATEDANNEALINGKERNEL_H_

#include "algorithms/SaUpdate.h"
#include "../util/rng/PhiloxWrapper.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
namespace culgt
{

namespace FullGaugeFixingSimulatedAnnealingKernel
{
	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelSaStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE temperature, int rngSeed, int rngCounter )
	{

#ifdef DPUPDATE
		typedef SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE, double> Algorithm;
#else
		typedef SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
#endif

		PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
		Algorithm saupdate( temperature, rng );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}

}


#endif
