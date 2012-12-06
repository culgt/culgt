/*
 * kernel_launch_bounds.h
 *
 *  Created on: Dez 06, 2012
 *      Author: schroeck
 * 
 * Here we define the values of the minBlocksPerMultiprocessor
 * in __launch_bounds__(256, minBlocksPerMultiprocessor)
 * for the kernels:
 * 
 * OR = overrelaxation
 * SR = stochastic relaxation
 * SA = simulated annealing
 * MS = micro step
 * 
 */

#ifndef KERNEL_LAUNCH_BOUNDS_H_
#define KERNEL_LAUNCH_BOUNDS_H_


#ifdef DOUBLEPRECISION
	#define OR_MINBLOCKS 1
	#define MS_MINBLOCKS 1
	#define SA_MINBLOCKS 1
	#define SR_MINBLOCKS 1
#else
	#define OR_MINBLOCKS 4
	#define MS_MINBLOCKS 4
	#define SA_MINBLOCKS 4
	#define SR_MINBLOCKS 4
#endif
	
	
#endif /* KERNEL_LAUNCH_BOUNDS_H_ */
