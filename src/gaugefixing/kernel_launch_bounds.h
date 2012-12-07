/*
 * kernel_launch_bounds.h
 *
 *  Created on: Dez 06, 2012
 *      Author: schroeck
 * 
 * Here we define the values of the minBlocksPerMultiprocessor
 * in __launch_bounds__(8*NSB, minBlocksPerMultiprocessor)
 * for the kernels:
 * 
 * OR = overrelaxation
 * SR = stochastic relaxation
 * SA = simulated annealing
 * MS = micro step
 * 
 *
 * 32 registers is optimal for OR.
 * with NSB = 32 this results in launch_bounds( threadsPerBlock=8*NSB=256, minBlocksPerSM=4 )
 * to keep register usage constant while changing threadsPerBlock, we have to adjust minBlocksPerSM
 *
 * Test launchbounds for DP updates in SP code (mixed code)
 */

#ifndef KERNEL_LAUNCH_BOUNDS_H_
#define KERNEL_LAUNCH_BOUNDS_H_

#include "GlobalConstants.h"

#ifdef DOUBLEPRECISION
	const int OR_MINBLOCKS = 1;
	const int MS_MINBLOCKS = 1;
	const int SA_MINBLOCKS = 1;
	const int SR_MINBLOCKS = 1;
#else

#ifdef USE_DP_ORUPDATE
	const int OR_MINBLOCKS = 128/NSB;
#else
	const int OR_MINBLOCKS = 128/NSB;
#endif

#ifdef USE_DP_SAUPDATE
	const int SA_MINBLOCKS = 128/NSB;
#else
	const int SA_MINBLOCKS = 128/NSB;
#endif

#ifdef USE_DP_MICROUPDATE
	const int MS_MINBLOCKS = 128/NSB;
#else
	const int MS_MINBLOCKS = 128/NSB;
#endif

#ifdef USE_DP_SRUPDATE
	const int SR_MINBLOCKS = 128/NSB;
#else
	const int SR_MINBLOCKS = 128/NSB;
#endif
	
#endif
	
#endif /* KERNEL_LAUNCH_BOUNDS_H_ */
