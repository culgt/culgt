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
 * This class contains the kernels and wrappers called
 * by MultiGPU_MPI_Communicator.hxx.
 * 
 * kernels as class members are not supported (even static): 
 * wrap the kernel calls and hide the kernels in namespace.
 * 
 */

#ifndef MULTIGPU_MPI_LANDAUKERNELSSU3_H_
#define MULTIGPU_MPI_LANDAUKERNELSSU3_H_

#include "./MultiGPU_MPI_AlgorithmOptions.h"


// kernels: 
namespace MPILKSU3
{
static const int Ndim = 4;
static const int Nc = 3;
template<class Algorithm> inline __device__ void applyOneTimeslice( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, Algorithm algorithm  );
__global__ void generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA );
__global__ void randomTrafo( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, int rngSeed, int rngCounter );
__global__ void orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter );
__global__ void microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity );
__global__ void saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter );
__global__ void projectSU3( Real* Ut );
__global__ void setHot( Real* Ut, int rngSeed, int rngCounter );
}

// wrappers:
class MultiGPU_MPI_LandauKernelsSU3
{
public:
	// constructor
	MultiGPU_MPI_LandauKernelsSU3();
	// tell CUDA to prefer the L1 cache
	static void initCacheConfig();
	// applies an anlgorithm (given in algoOptions) to a single timeslice
	void applyOneTimeslice( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, MultiGPU_MPI_AlgorithmOptions algoOptions  );
	// projects all SU(3) matrices in a timeslice back to the group
	void projectSU3( int a, int b, cudaStream_t stream, Real* Ut );
	// fill gauge field on the devices with random SU(3) matrices
	void setHot( int a, int b, cudaStream_t stream, Real* Ut, MultiGPU_MPI_AlgorithmOptions algoOptions );
	// generates the gauge quality on a timesclice for all sites, no reduction
	void generateGaugeQualityPerSite( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA );
	
private:
	
};



#endif /* MULTIGPU_MPI_LANDAUKERNELSSU3_H_ */
