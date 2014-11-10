/*
 * LaunchBounds.h
 *
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef LAUNCHBOUNDS_H_
#define LAUNCHBOUNDS_H_

template<int TmaxThreadsPerBlock, int TminBlocksPerMultiprocessor> struct LaunchBounds
{
	static const int maxThreadsPerBlock = TmaxThreadsPerBlock;
	static const int minBlocksPerMultiprocessor = TminBlocksPerMultiprocessor;
};


#endif /* LAUNCHBOUNDS_H_ */
