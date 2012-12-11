/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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
