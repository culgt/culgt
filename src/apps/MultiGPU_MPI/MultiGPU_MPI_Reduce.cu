/**
 * Reduce.hxx
 *
 *  Created on: Nov. 30, 2012
 *      Author: schroeck
 *
 */

#include "./MultiGPU_MPI_Reduce.h"
#include "../../util/cuda/cuda_host_device.h"

/**
 * TODO place somewhere else
 */
CUDA_HOST_DEVICE bool isPowerOfTwo( unsigned long x )
{
    return (x & (x - 1)) == 0;
}

// __device__ inline double cuFabs( double a )
// {
// 	return (a>0)?(a):(-a);
// }




Reduce::Reduce( int _size )
{
	size = _size;
	initReductionBlockSize();
}

Reduce::~Reduce()
{
}

/**
 * Choose red. block size such that latticesize/redBlockSize/redBlockSize is an integer and redBlockSize is a power of 2
 */
void Reduce::initReductionBlockSize()
{
	int pot = 1;
	for( int i = 1; i < 10; i++ )
	{
		pot *= 2;
		if( size % (pot*pot) == 0 )
		{
			redBlockSize = pot;
		}
		else
		{
			break;
		}
	}
// 	printf( "We chose reduction block size = %d\n", redBlockSize );
// 	if( isPowerOfTwo( size/redBlockSize/redBlockSize) ) printf( "We use parallel reduction in last reduction step\n" );
// 	else printf( "We can't use parallel reduction in last reduction step\n" );
}

void Reduce::reduceStep1( int a, int b, cudaStream_t stream, double *Ar )
{
	MPIRED::reduceStep1<<<a,b,0,stream>>>( Ar );
}

void Reduce::reduceStep2( int a, int b, cudaStream_t stream, double *Ar )
{
	MPIRED::reduceStep2<<<a,b,0,stream>>>( Ar );
}

void Reduce::reduceStep3( int a, int b, cudaStream_t stream, double *Ar, int redBlockSize )
{
	MPIRED::reduceStep3<<<a,b,0,stream>>>( Ar, redBlockSize );
}

void Reduce::reduceSerial( int a, int b, cudaStream_t stream, double *Ar )
{
	MPIRED::reduceSerial<<<a,b,0,stream>>>( Ar, size );
}

double Reduce::getReducedValue( cudaStream_t stream, double *Ar )
{
	MPIRED::reduceStep1<<<size/redBlockSize,redBlockSize,0,stream>>>( Ar );
	MPIRED::reduceStep2<<<size/redBlockSize/redBlockSize,redBlockSize,0,stream>>>( Ar );
	MPIRED::reduceStep3<<<1,size/redBlockSize/redBlockSize,0,stream>>>( Ar, redBlockSize );

// 	MPIRED::reduceSerial<<<1,1,0,stream>>>( Ar, size );

	cudaMemcpy( &currentAr, Ar, sizeof(double), cudaMemcpyDeviceToHost );
	
	return currentAr;
}

void Reduce::setArrayZero( cudaStream_t stream, double *Ar )
{
	MPIRED::setArrayZero<<<size/redBlockSize,redBlockSize,0,stream>>>( Ar );
}

// kernels: hide in namespace
namespace MPIRED
{
__global__ void reduceStep1( double *Ar )
{
	int site = blockIdx.x*blockDim.x + threadIdx.x;
	
	// reduce within each thread block
	for( int stride=blockDim.x/2; stride>0; stride>>=1 )
	{
		__syncthreads();
		if( threadIdx.x<stride ) 
		{
			Ar[site] += Ar[site+stride];
		}
	}
}

__global__ void reduceStep2( double *Ar )
{
	int i = blockIdx.x*blockDim.x*blockDim.x + threadIdx.x*blockDim.x;
	
	for( int stride=blockDim.x/2; stride>0; stride>>=1 )
	{
		__syncthreads();
		if( threadIdx.x<stride ) 
		{
			Ar[i] += Ar[i+stride*blockDim.x];
		}
	}
}


__global__ void reduceStep3( double *Ar, int redBlockSize )
{
	int bs = redBlockSize*redBlockSize;
	int i = bs*threadIdx.x;

	if( isPowerOfTwo( blockDim.x ) ) // we can use parallel reduction
	{
		for( int stride=blockDim.x/2; stride>0; stride>>=1 )
		{
			__syncthreads();
			if( threadIdx.x<stride )
			{
				Ar[i] += Ar[i+stride*bs];
			}
		}
	}
	else // serialise the last iteration (only one thread is active)
	{
		if (i == 0 )
			for( int j = 1; j < blockDim.x; j ++ )
			{
				Ar[0] += Ar[0+bs*j];
			}
	}
}


__global__ void reduceSerial( double *Ar, int size )
{
	double A = 0;
	for( int i = 0; i < size; i++ )
	{
		A+= Ar[i];
	}
	Ar[0] = A;
}

__global__ void setArrayZero( double *Ar )
{
	Ar[blockIdx.x*blockDim.x + threadIdx.x] = 0.0;
}

// __device__ double average_or_max( double a, double b, StoppingCrit ma)
// {
// 	if( ma == AVERAGE ) return a+b;
// 	else if ( ma == MAX ) return (a>b?a:b);
// 	else return 0;
// }
}

