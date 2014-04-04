/**
 * SubgroupIterator.h
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef SUBGROUPITERATOR_H_
#define SUBGROUPITERATOR_H_

#include "../cudacommon/cuda_host_device.h"

template<int Nc> class SubgroupIterator
{
};

template<> class SubgroupIterator<3>
{
public:
	static const int NSubgroups = 3;
	template<typename SubgroupStep> CUDA_HOST_DEVICE static inline void iterate( SubgroupStep& aSubgroupStep )
	{
		// the order of subgroups is important for performance! try all combinations, i.e. permutations of permutations of (0,1),(0,2),(1,2)
		aSubgroupStep.subgroupStep( 0, 2 );
		aSubgroupStep.subgroupStep( 1, 2 );
		aSubgroupStep.subgroupStep( 1, 0 );
	}
};

template<> class SubgroupIterator<2>
{
public:
	static const int NSubgroups = 1;
	template<typename SubgroupStep> CUDA_HOST_DEVICE static inline void iterate( SubgroupStep& aSubgroupStep )
	{
		aSubgroupStep.subgroupStep( 0, 1 );
	}
};


#endif /* SUBGROUPITERATOR_H_ */
