/**
 * GaugeConfigurationCudaHelper.h
 *
 *  Created on: Feb 26, 2014
 *      Author: vogt
 */

#ifndef GAUGECONFIGURATIONCUDAHELPER_H_
#define GAUGECONFIGURATIONCUDAHELPER_H_

#include "../cudacommon/cuda_error.h"


template<typename T> class GaugeConfigurationCudaHelper
{
public:
	static void allocateMemory( T** pointer, size_t size );
	static void freeMemory( T* pointer );
	static void setElement( T* pointer, int i, const T val );
	template<typename ConfigurationPattern, typename LinkType> static void setLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu, const LinkType link  );
	template<typename ConfigurationPattern, typename LinkType> static LinkType getLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu );
	static T getElement( T* pointer, int i );
	static void copyToDevice( T* dest, T* src, size_t size );
	static void copyToHost( T* dest, T* src, size_t size );
};

template<typename T> void GaugeConfigurationCudaHelper<T>::allocateMemory( T** pointerToPointer, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, configurationSize*sizeof(T) ) , "cudaMalloc in GaugeConfigurationHelper" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::freeMemory( T* pointer )
{
	CUDA_SAFE_CALL( cudaFree( pointer ), "cudaFree in GaugeConfigurationHelper" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::setElement( T* pointer, int i, T val )
{
	CUDA_SAFE_CALL( cudaMemcpy( &pointer[i], &val, sizeof(T), cudaMemcpyHostToDevice ) , "cudaMemcpy in GaugeConfigurationHelper<T>::setElement" );
}

template<typename T> template<typename ConfigurationPattern, typename LinkType> void GaugeConfigurationCudaHelper<T>::setLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu, const LinkType link  )
{
	for( int i = 0; i < ConfigurationPattern::PARAMTYPE::SIZE; i++ )
	{
		GaugeConfigurationCudaHelper<T>::setElement( pointer, ConfigurationPattern::getIndex(site,mu,i), link.get(i) );
	}
}

template<typename T> T GaugeConfigurationCudaHelper<T>::getElement( T* pointer, int i )
{
	T val;
	CUDA_SAFE_CALL( cudaMemcpy( &val, &pointer[i], sizeof(T), cudaMemcpyDeviceToHost ) , "cudaMemcpy in GaugeConfigurationHelper<T>::getElement" );
	return val;
}

template<typename T> template<typename ConfigurationPattern, typename LinkType> LinkType GaugeConfigurationCudaHelper<T>::getLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu )
{
	LinkType link;
	for( int i = 0; i < ConfigurationPattern::PARAMTYPE::SIZE; i++ )
	{
		T element = GaugeConfigurationCudaHelper<T>::getElement( pointer, ConfigurationPattern::getIndex(site,mu,i) );
		link.set( i, element );
	}
	return link;
}

template<typename T> void GaugeConfigurationCudaHelper<T>::copyToDevice( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyHostToDevice ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToDevice" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::copyToHost( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyDeviceToHost ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToHost" );
}


#endif /* GAUGECONFIGURATIONHELPER_H_ */
