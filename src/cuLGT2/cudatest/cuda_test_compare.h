/**
 * cuda_test_compare.h
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef CUDA_TEST_COMPARE_H_
#define CUDA_TEST_COMPARE_H_

#include "gmock/gmock.h"
#include "../../lattice/cuda/CudaError.h"

template<typename Kernel, typename InputType, typename OutputType> __global__ void cudaMetaKernel( const Kernel kernel, const InputType& in, OutputType& out )
{
	kernel( in, out );
}

template<typename Kernel, typename OutputType> __global__ void cudaMetaKernel( const Kernel kernel, OutputType& out )
{
	kernel( out );
}


template<typename Kernel, typename InputType, typename OutputType>  void runOnDevice( const Kernel kernel, const InputType& in, OutputType& out )
{
	InputType* dIn;
	cudaMalloc( (void**)(&dIn), sizeof( InputType ) );

	OutputType* dOut;
	cudaMalloc( (void**)(&dOut), sizeof( OutputType ) );
	CudaError::getLastError( "cudaMalloc" );

	cudaMemcpy( dIn, &in, sizeof( InputType ), cudaMemcpyHostToDevice );
	CudaError::getLastError( "memcopy to device" );


	cudaMetaKernel<<<1,1>>>(kernel, *dIn, *dOut);
	CudaError::getLastError( "kernel call" );


	cudaMemcpy( &out, dOut, sizeof( OutputType ), cudaMemcpyDeviceToHost );
	CudaError::getLastError( "memcopy to host" );
}

template<typename Kernel, typename OutputType>  void runOnDevice( const Kernel kernel,OutputType& out )
{
	OutputType* dOut;
	cudaMalloc( (void**)(&dOut), sizeof( OutputType ) );
	CudaError::getLastError( "cudaMalloc" );


	cudaMetaKernel<<<1,1>>>( kernel, *dOut);
	CudaError::getLastError( "kernel call" );


	cudaMemcpy( &out, dOut, sizeof( OutputType ), cudaMemcpyDeviceToHost );
	CudaError::getLastError( "memcopy to host" );
}

template<typename Kernel, typename InputType, typename OutputType>  void runOnCpu( const Kernel kernel, const InputType& in, OutputType& out )
{
	kernel( in, out );
}

template<typename Kernel, typename OutputType>  void runOnCpu( const Kernel kernel, OutputType& out )
{
	kernel( out );
}

template<typename Kernel, typename InputType, typename OutputType> void cudaRunAndCompare( const Kernel& kernel, InputType& in, OutputType& out, bool expect = true )
{
	runOnCpu( kernel, in, out );

	OutputType cudaOut;
	runOnDevice( kernel, in, cudaOut );

	ASSERT_THAT( cudaOut.isEqual( out ), expect );
}

template<typename Kernel, typename OutputType> void cudaRunAndCompare( const Kernel& kernel, OutputType& out, bool expect = true )
{
	runOnCpu( kernel, out );

	OutputType cudaOut;
	runOnDevice( kernel, cudaOut );

	ASSERT_THAT( cudaOut.isEqual( out ), expect );
}

#endif /* CUDA_TEST_COMPARE_H_ */
