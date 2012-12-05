/**
 * Reduce.hxx
 *
 *  Created on: Nov. 30, 2012
 *      Author: schroeck
 * 
 * This class and the kernels in the namespace
 * reduces (sums over) an array Ar of doubles to a single double
 * by calling three levels of tree-like reducers.
 *
 */

#ifndef MULTIGPU_MPI_REDUCE_HXX_
#define MULTIGPU_MPI_REDUCE_HXX_


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
