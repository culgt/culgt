/**
 * RunInfo.h
 *
 *  Created on: Mar 26, 2014
 *      Author: vogt
 */

#ifndef RUNINFO_H_
#define RUNINFO_H_

#include "../lattice/SubgroupIterator.h"

class RunInfo
{
public:
	RunInfo(): gflops(0), throughput(0)
	{
	}
	RunInfo( double gflops, double throughput ): gflops(gflops), throughput(throughput)
	{
	}

	double getGflops() const
	{
		return gflops;
	}

	double getThroughput() const
	{
		return throughput;
	}


	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType> static RunInfo makeRunInfo( int sites, double time, int iter, int flopsForAlgorithm )
	{
		int NSumOfLocalUpdate = GaugeType::LinksInvolved * ( LocalLinkType::PARAMTYPE::FlopsGetSU2Subgroup + 4 ); // #links * ( ... + additions in gather )
		int NApplyTrafo = 8*LocalLinkType::PARAMTYPE::FlopsSubgroupMult; // #links * subgroupMult

		long flopsPerSite = SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::NSubgroups * ( NSumOfLocalUpdate + flopsForAlgorithm + NApplyTrafo);// #subgroup*( sumOfLocalUpdate + updatealgorithm + applyTrafo )

		long bytesPerSite = 2*8*GlobalLinkType::PATTERNTYPE::PARAMTYPE::SIZE*sizeof(typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE); // read/write * #links(NDim*2) * sizeof(link)

		double gflops = (double)(flopsPerSite*(long)sites*(long)iter)/time/1.0e9;
		double throughput = (double)(bytesPerSite*(long)sites*(long)(iter))/time/1.0e9;

		return RunInfo( gflops, throughput );
	}

private:
	double gflops;
	double throughput;
};


#endif /* RUNINFO_H_ */
