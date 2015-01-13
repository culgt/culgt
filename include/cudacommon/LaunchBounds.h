/*
 * LaunchBounds.h
 *
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef LAUNCHBOUNDS_H_
#define LAUNCHBOUNDS_H_

#include <iostream>
#include <string>
using std::ostream;
using std::string;

template<int TmaxThreadsPerBlock, int TminBlocksPerMultiprocessor> struct LaunchBounds
{
	static const int maxThreadsPerBlock = TmaxThreadsPerBlock;
	static const int minBlocksPerMultiprocessor = TminBlocksPerMultiprocessor;

	static LaunchBounds<TmaxThreadsPerBlock,TminBlocksPerMultiprocessor> printable()
	{
		return *(new LaunchBounds<TmaxThreadsPerBlock,TminBlocksPerMultiprocessor>);
	}
};

template<int TmaxThreadsPerBlock, int TminBlocksPerMultiprocessor> ostream& operator<<(ostream& out, const LaunchBounds<TmaxThreadsPerBlock,TminBlocksPerMultiprocessor>& lb)
{
	out << "LB(" << lb.maxThreadsPerBlock << "," << lb.minBlocksPerMultiprocessor << ")";
    return out;
}

#endif /* LAUNCHBOUNDS_H_ */
