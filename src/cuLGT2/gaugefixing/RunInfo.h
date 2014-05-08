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

	void print()
	{
		std::cout << std::fixed << std::setprecision( 2 ) << gflops << " GFlop/s at " << throughput << " GBytes/s memory throughput." << std::endl;
	}

	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType> static RunInfo makeRunInfo( int sites, double time, long iter, int flopsForAlgorithm )
	{
		int NSumOfLocalUpdate = GaugeType::LinksInvolved * ( LocalLinkType::PARAMTYPE::FlopsGetSU2Subgroup + 4 ); // #links * ( ... + additions in gather )
		int NApplyTrafo = 8*LocalLinkType::PARAMTYPE::FlopsSubgroupMult; // #links * subgroupMult

		long flopsPerSite = SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::NSubgroups * ( NSumOfLocalUpdate + flopsForAlgorithm + NApplyTrafo);// #subgroup*( sumOfLocalUpdate + updatealgorithm + applyTrafo )

		long bytesPerSite = 2*8*GlobalLinkType::PATTERNTYPE::PARAMTYPE::SIZE*sizeof(typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE); // read/write * #links(NDim*2) * sizeof(link)

		double gflops = (double)(flopsPerSite*(long)sites*(long)iter)/time/1.0e9;
		double throughput = (double)(bytesPerSite*(long)sites*(long)(iter))/time/1.0e9;

		return RunInfo( gflops, throughput );
	}

	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType> static RunInfo makeRunInfo( int sites, double time, long iter1, int flopsForAlgorithm1, long iter2, int flopsForAlgorithm2 )
	{
		int NSumOfLocalUpdate = GaugeType::LinksInvolved * ( LocalLinkType::PARAMTYPE::FlopsGetSU2Subgroup + 4 ); // #links * ( ... + additions in gather )
		int NApplyTrafo = 8*LocalLinkType::PARAMTYPE::FlopsSubgroupMult; // #links * subgroupMult

		long flopsPerSite1 = SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::NSubgroups * ( NSumOfLocalUpdate + flopsForAlgorithm1 + NApplyTrafo);// #subgroup*( sumOfLocalUpdate + updatealgorithm + applyTrafo )
		long flopsPerSite2 = SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::NSubgroups * ( NSumOfLocalUpdate + flopsForAlgorithm2 + NApplyTrafo);// #subgroup*( sumOfLocalUpdate + updatealgorithm + applyTrafo )

		long bytesPerSite = 2*8*GlobalLinkType::PATTERNTYPE::PARAMTYPE::SIZE*sizeof(typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE); // read/write * #links(NDim*2) * sizeof(link)

		double gflops = (double)((flopsPerSite1*iter1+flopsPerSite2*iter2)*(long)sites)/time/1.0e9;
		double throughput = (double)(bytesPerSite*(long)sites*(iter1+iter2))/time/1.0e9;

		return RunInfo( gflops, throughput );
	}

private:
	double gflops;
	double throughput;
};


#endif /* RUNINFO_H_ */
