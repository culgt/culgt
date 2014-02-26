/**
 * GaugeConfiguration.cu
 *
 * @note I separated implementation and header of this template because nvcc has problems with exceptions. That's why I do not want any CUDA code in GaugeConfiguration.
 * 		The explicit template instantiations are a big drawback and a potential source of error.
 *
 *  Created on: Feb 26, 2014
 *      Author: vogt
 */

#include "../cudacommon/cuda_error.h"
#include "GaugeConfigurationHelper.h"

/**
 *
 * @param pointer
 * @param configurationSize number of T-elements in the array (NOT number of bytes!)
 */
template<typename T> void GaugeConfigurationHelper<T>::allocateMemory( T** pointerToPointer, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, configurationSize*sizeof(T) ) , "cudaMalloc in GaugeConfigurationHelper" );
}

template<typename T> void GaugeConfigurationHelper<T>::freeMemory( T* pointer )
{
	CUDA_SAFE_CALL( cudaFree( pointer ), "cudaFree in GaugeConfigurationHelper" );
}

template<typename T> void GaugeConfigurationHelper<T>::setElement( T* pointer, int i, T val )
{
	CUDA_SAFE_CALL( cudaMemcpy( &pointer[i], &val, sizeof(T), cudaMemcpyHostToDevice ) , "cudaMemcpy in GaugeConfigurationHelper<T>::setElement" );
}

template<typename T> T GaugeConfigurationHelper<T>::getElement( T* pointer, int i )
{
	T val;
	CUDA_SAFE_CALL( cudaMemcpy( &val, &pointer[i], sizeof(T), cudaMemcpyDeviceToHost ) , "cudaMemcpy in GaugeConfigurationHelper<T>::getElement" );
	return val;
}

template<typename T> void GaugeConfigurationHelper<T>::copyToDevice( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyHostToDevice ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToDevice" );
}

template<typename T> void GaugeConfigurationHelper<T>::copyToHost( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyDeviceToHost ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToHost" );
}

template class GaugeConfigurationHelper<float>;
template class GaugeConfigurationHelper<double>;
