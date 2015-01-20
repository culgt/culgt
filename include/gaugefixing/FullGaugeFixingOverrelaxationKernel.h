/*
 * FullGaugeFixingOverrelaxationKernel.h
 *
 *  Created on: Nov 4, 2014
 *      Author: vogt
 */

#ifndef FULLGAUGEFIXINGOVERRELAXATIONKERNEL_H_
#define FULLGAUGEFIXINGOVERRELAXATIONKERNEL_H_

#include "algorithms/MicroUpdate.h"
#include "algorithms/OrUpdate.h"

namespace culgt
{


namespace FullGaugeFixingOverrelaxationKernel
{
	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
	{
#ifdef DPUPDATE
		typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE, double> Algorithm;
#else
		typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
#endif

		Algorithm orupdate( orParameter );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}

	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelMicroStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity )
	{
#ifdef DPUPDATE
		typedef MicroUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE, double> Algorithm;
#else
		typedef MicroUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
#endif

		Algorithm microupdate;

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( microupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( microupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}

}


#endif /* FULLGAUGEFIXINGOVERRELAXATIONKERNEL_H_ */
