/*
 * TimesliceGaugeFixingOverrelaxationKernel.h
 *
 *  Created on: Nov 4, 2014
 *      Author: vogt
 */

#ifndef TIMESLICEGAUGEFIXINGSIMULATEDANNEALINGKERNEL_H_
#define TIMESLICEGAUGEFIXINGSIMULATEDANNEALINGKERNEL_H_
#include "algorithms/SaUpdate.h"

namespace culgt
{

namespace TimesliceGaugeFixingSimulatedAnnealingKernel
{
	template<typename GlobalLinkType,typename GlobalLinkType2, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelSaStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE temperature, int rngSeed, int rngCounter )
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
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}

}


#endif
