/**
 * Reduce.hxx
 *
 *  Created on: Nov. 30, 2012
 *      Author: schroeck
 *
 * TODO can this go in one file .hxx ?
 */

#ifndef MULTIGPU_MPI_REDUCE_HXX_
#define MULTIGPU_MPI_REDUCE_HXX_

// #include "GlobalConstants.hxx"
// #include "../../util/datatype/datatypes.h"
// #include "../../util/cuda/cuda_host_device.h"
// #include "../access_pattern/StandardPattern.hxx"
// #include "../access_pattern/GpuCoulombPattern.hxx"
// #include "../access_pattern/GpuLandauPattern.hxx"
// #include "../SiteCoord.hxx"
// #include "../SiteIndex.hxx"
// #include "../Link.hxx"
// #include "../SU3.hxx"
// #include "../SU2.hxx"
// #include "../Matrix.hxx"
// #include "../LinkFile.hxx"

// #include <iostream>



// kernels: hide in namespace
namespace MPIRED
{
__global__ void reduceStep1( double *Ar );
__global__ void reduceStep2( double *Ar );
__global__ void reduceStep3( double *Ar, int redBlockSize );
__global__ void reduceSerial( double *Ar, int size );
__global__ void setArrayZero( double *Ar );
// __device__ double average_or_max( double a, double b );
}

class Reduce
{
public:
	Reduce( int _size );
	~Reduce();
	double getReducedValue( cudaStream_t stream, double *Ar );
	void setArrayZero( cudaStream_t stream, double *Ar );
private:
	int size;
	int redBlockSize;
	double currentAr;
	void initReductionBlockSize();
	// kernel wrapper
	void reduceStep1( int a, int b, cudaStream_t stream, double *Ar );
	void reduceStep2( int a, int b, cudaStream_t stream, double *Ar );
	void reduceStep3( int a, int b, cudaStream_t stream, double *Ar, int redBlockSize );
	void reduceSerial( int a, int b, cudaStream_t stream, double *Ar );
};


#endif /* MULTIGPU_MPI_REDUCE_HXX_ */
