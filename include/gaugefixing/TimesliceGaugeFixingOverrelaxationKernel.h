/*
 * TimesliceGaugeFixingOverrelaxationKernel.h
 *
 *  Created on: Nov 4, 2014
 *      Author: vogt
 */

#ifndef TIMESLICEGAUGEFIXINGOVERRELAXATIONKERNEL_H_
#define TIMESLICEGAUGEFIXINGOVERRELAXATIONKERNEL_H_

#include "algorithms/MicroUpdate.h"
#include "algorithms/OrUpdate.h"

namespace culgt
{

namespace TimesliceGaugeFixingOverrelaxationKernel
{
	template<typename GlobalLinkType,typename GlobalLinkType2, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter)
	{
		typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
		Algorithm orupdate( orParameter );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}

	template<typename GlobalLinkType,typename GlobalLinkType2, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelMicroStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity )
	{
		typedef MicroUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
		Algorithm microupdate;

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( microupdate );
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( microupdate );
			gaugefixing.template applyAlgorithmTimeslice<GlobalLinkType2>( Ut, UtDown, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}

}


#endif /* TIMESLICEGAUGEFIXINGOVERRELAXATIONKERNEL_H_ */
