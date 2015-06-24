/**
 * GaugeConfigurationCudaHelper.h
 *
 *  Created on: Feb 26, 2014
 *      Author: vogt
 */

#ifndef GAUGECONFIGURATIONCUDAHELPER_H_
#define GAUGECONFIGURATIONCUDAHELPER_H_

#include "lattice/site_indexing/ParityType.h"
#include "cudacommon/cuda_error.h"
#include "LatticeDimension.h"
#include "LocalLink.h"
#include "GlobalLink.h"
#include "parameterization_types/SU2Vector4.h"
#include "parameterization_types/SUNRealFull.h"
#include "parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"
#include "parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "SubgroupIterator.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "parameterization_types/SUNRealFull.h"
#include "KernelSetup.h"
#include "cudacommon/DeviceMemoryManager.h"
#include <boost/mpl/assert.hpp>

namespace culgt
{

template<typename T> class GaugeConfigurationCudaHelper
{
public:
	static void allocateMemory( T** pointer, size_t size );
	static void freeMemory( T* pointer );
	static void setElement( T* pointer, int i, const T val );
	template<typename ConfigurationPattern, typename LinkType> static void setLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu, const LinkType link  );
	template<typename ConfigurationPattern, typename LinkType> static LinkType getLink( T* pointer, const typename ConfigurationPattern::SITETYPE site, const int mu );
	template<typename ConfigurationPattern, typename RNG> static void setHot( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter );
	template<typename ConfigurationPattern, typename RNG> static void setRandomZN( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter, float percentage );
	template<typename ConfigurationPattern, typename LocalLinkType, typename RNG> static void reproject( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter );
	template<typename ConfigurationPattern> static void setCold( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim );
	static T getElement( T* pointer, int i );
	static void copyDeviceToDevice( T* dest, T* src, size_t size );
	static void copyToDevice( T* dest, T* src, size_t size );
	static void copyToHost( T* dest, T* src, size_t size );
};

template<typename T> void GaugeConfigurationCudaHelper<T>::allocateMemory( T** pointerToPointer, size_t configurationSize )
{
	culgt::DeviceMemoryManager::malloc( pointerToPointer, configurationSize*sizeof(T), "Gauge configuration" );
//	CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, configurationSize*sizeof(T) ) , "cudaMalloc in GaugeConfigurationHelper" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::freeMemory( T* pointer )
{
	culgt::DeviceMemoryManager::free( pointer );
//	CUDA_SAFE_CALL( cudaFree( pointer ), "cudaFree in GaugeConfigurationHelper" );
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


template<typename T> void GaugeConfigurationCudaHelper<T>::copyDeviceToDevice( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyDeviceToDevice ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyDeviceToDevice" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::copyToDevice( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyHostToDevice ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToDevice" );
}

template<typename T> void GaugeConfigurationCudaHelper<T>::copyToHost( T* dest, T* src, size_t configurationSize )
{
	CUDA_SAFE_CALL( cudaMemcpy( dest, src, configurationSize*sizeof(T), cudaMemcpyDeviceToHost ), "cudaMemcpy in GaugeConfigurationHelper<T>::copyToHost" );
}

#ifdef __CUDACC__

namespace GaugeConfigurationCudaHelperKernel
{
	template<typename LocalLinkType, typename RNG> class SetHotSubgroup
	{
	public:
		__device__ SetHotSubgroup( LocalLinkType& link, RNG& rng ): link(link), rng(rng) {};
		__device__ void subgroupStep( int i, int j )
		{
			culgt::LocalLink<culgt::SU2Vector4<typename LocalLinkType::PARAMTYPE::REALTYPE> > quaternion = link.getSU2Subgroup( i, j );
			typename Real4<typename LocalLinkType::PARAMTYPE::REALTYPE>::VECTORTYPE& q = quaternion[0];

			typename LocalLinkType::PARAMTYPE::REALTYPE alpha, phi, cos_theta, sin_theta, sin_alpha;
			alpha = rng.rand();
			phi = 2.0 * rng.rand();
			cos_theta = 2.0 * rng.rand() - 1.0;
			sin_theta = ::sqrt(1.0 - cos_theta * cos_theta);
			sin_alpha = sinpi(alpha);
			q.x = cospi(alpha);
			q.y = sin_alpha * sin_theta * cospi(phi);
			q.z = sin_alpha * sin_theta * sinpi(phi);
			q.w = sin_alpha * cos_theta;

			quaternion[0] = q;
			link.rightSubgroupMult( quaternion, i, j );
		}
	private:
		LocalLinkType& link;
		RNG& rng;
	};


	template<typename ConfigurationPattern, typename LocalLinkType, typename RNG> __global__ void reproject( typename ConfigurationPattern::PARAMTYPE::TYPE* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter )
	{
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		RNG rng( index, rngSeed, rngCounter );

		typename ConfigurationPattern::SITETYPE site( dim, DO_NOT_USE_NEIGHBOURS );
		site.setIndex( index );
		for( int mu = 0; mu < ConfigurationPattern::SITETYPE::NDIM; mu++ )
		{
			LocalLinkType local;
			culgt::GlobalLink<ConfigurationPattern> global( pointer, site, mu );

			local = global;

			local.reproject();

			global = local;
		}
	}

	/**
	 * FIXME: why not use the same template parameters as reproject??
	 * @param pointer
	 * @param dim
	 * @param rngSeed
	 * @param rngCounter
	 */
	template<typename T, typename ConfigurationPattern, typename RNG> __global__ void setHot( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter )
	{
		typedef culgt::LocalLink<culgt::SUNRealFull<ConfigurationPattern::PARAMTYPE::NC, typename ConfigurationPattern::PARAMTYPE::REALTYPE> > LocalLinkType;

		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		RNG rng( index, rngSeed, rngCounter );

		typename ConfigurationPattern::SITETYPE site( dim, DO_NOT_USE_NEIGHBOURS );
		site.setIndex( index );
		for( int mu = 0; mu < ConfigurationPattern::SITETYPE::NDIM; mu++ )
		{
			LocalLinkType local;
			local.identity();

			SetHotSubgroup<LocalLinkType,RNG> aSetHotSubgroup( local, rng );

			SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( aSetHotSubgroup );

			culgt::GlobalLink<ConfigurationPattern> global( pointer, site, mu );

			global = local;
		}
	}

	/**
	 *
	 */
	template<typename T, typename ConfigurationPattern, typename RNG> __global__ void setRandomZN( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter, float percentage )
	{
		// TODO this is only Z(2)
		BOOST_MPL_ASSERT_RELATION( ConfigurationPattern::PARAMTYPE::NC, ==, 2 );

		typedef culgt::LocalLink<culgt::SUNRealFull<ConfigurationPattern::PARAMTYPE::NC, typename ConfigurationPattern::PARAMTYPE::REALTYPE> > LocalLinkType;

		int index = blockIdx.x * blockDim.x + threadIdx.x;
		VERIFY_LATTICE_SIZE( dim, index );

		RNG rng( index, rngSeed, rngCounter );

		typename ConfigurationPattern::SITETYPE site( dim, DO_NOT_USE_NEIGHBOURS );
		site.setIndex( index );
		for( int mu = 0; mu < ConfigurationPattern::SITETYPE::NDIM; mu++ )
		{
			LocalLinkType local;
			local.identity();
			if( rng.rand() < percentage )
				local *= -1.;

			culgt::GlobalLink<ConfigurationPattern> global( pointer, site, mu );
			global = local;
		}
	}

	template<typename T, typename ConfigurationPattern> __global__ void setCold( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim )
	{
		typedef culgt::LocalLink<culgt::SUNRealFull<ConfigurationPattern::PARAMTYPE::NC, typename ConfigurationPattern::PARAMTYPE::REALTYPE> > LocalLinkType;

		int index = blockIdx.x * blockDim.x + threadIdx.x;
		VERIFY_LATTICE_SIZE( dim, index );

		typename ConfigurationPattern::SITETYPE site( dim, DO_NOT_USE_NEIGHBOURS );
		site.setIndex( index );
		for( int mu = 0; mu < ConfigurationPattern::SITETYPE::NDIM; mu++ )
		{
			LocalLinkType local;
			local.identity();

			culgt::GlobalLink<ConfigurationPattern> global( pointer, site, mu );
			global = local;
		}
	}
}

template<typename T> template<typename ConfigurationPattern, typename LocalLinkType, typename RNG> void GaugeConfigurationCudaHelper<T>::reproject( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter )
{
	CUDA_LAST_ERROR( "before reproject" );
	KernelSetup<ConfigurationPattern::SITETYPE::NDIM> setup(dim);
	GaugeConfigurationCudaHelperKernel::reproject<ConfigurationPattern,LocalLinkType,RNG><<<setup.getGridSize(),setup.getBlockSize()>>>( pointer, dim, rngSeed, rngCounter );
	CUDA_LAST_ERROR( "kernel reproject" );
}

template<typename T> template<typename ConfigurationPattern, typename RNG> void GaugeConfigurationCudaHelper<T>::setHot( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter )
{
	KernelSetup<ConfigurationPattern::SITETYPE::NDIM> setup(dim);
	GaugeConfigurationCudaHelperKernel::setHot<T,ConfigurationPattern,RNG><<<setup.getGridSize(),setup.getBlockSize()>>>( pointer, dim, rngSeed, rngCounter );
	CUDA_LAST_ERROR( "kernel setHot" );
}

template<typename T> template<typename ConfigurationPattern, typename RNG> void GaugeConfigurationCudaHelper<T>::setRandomZN( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim, int rngSeed, int rngCounter, float percentage )
{
	KernelSetup<ConfigurationPattern::SITETYPE::NDIM> setup(dim);
	GaugeConfigurationCudaHelperKernel::setRandomZN<T,ConfigurationPattern,RNG><<<setup.getGridSize(),setup.getBlockSize()>>>( pointer, dim, rngSeed, rngCounter, percentage );
	CUDA_LAST_ERROR( "kernel setRandomZN" );
}

template<typename T> template<typename ConfigurationPattern> void GaugeConfigurationCudaHelper<T>::setCold( T* pointer, LatticeDimension<ConfigurationPattern::SITETYPE::NDIM> dim )
{
	KernelSetup<ConfigurationPattern::SITETYPE::NDIM> setup(dim);
	GaugeConfigurationCudaHelperKernel::setCold<T,ConfigurationPattern><<<setup.getGridSize(),setup.getBlockSize()>>>( pointer, dim );
	CUDA_LAST_ERROR( "kernel setCold" );
}

#endif

}

#endif /* GAUGECONFIGURATIONHELPER_H_ */
