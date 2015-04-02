/**
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef WILSONLOOPAVERAGECUDAHELPER_H_
#define WILSONLOOPAVERAGECUDAHELPER_H_

#include "lattice/LatticeDimension.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "cudacommon/cuda_error.h"
#include "lattice/KernelSetup.h"
#include "lattice/site_indexing/SiteNeighbourTableManager.h"
#include "WilsonLoop.h"

namespace culgt
{

namespace WilsonLoopAverageCudaHelperKernel
{
	template<typename PatternType, typename LocalLinkType> __global__ void kernelCalculateWilsonLoops( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn, int dir1, int length1, int dir2, int length2 )
	{
		typedef typename PatternType::PARAMTYPE::REALTYPE T;
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site(dim, nn );
		site.setLatticeIndex( index );

		WilsonLoop<PatternType,LocalLinkType> w( U, dim );

		data[index] = w.getReTraceWilsonLoopNormalized( site,  dir1, length1, dir2, length2 );

	};

	template<typename PatternType, typename LocalLinkType> __global__ void kernelCalculatePolyakovLoopCorrelators( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn, int dirSeparation, int separation, int timeDir = 0 )
	{
		typedef typename PatternType::PARAMTYPE::REALTYPE T;
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site(dim, nn );
		site.setLatticeIndex( index );

		WilsonLoop<PatternType,LocalLinkType> w( U, dim );

		data[index] = w.getPolyakovLoopCorrelator( site, dirSeparation, separation, timeDir );

	};

	template<typename PatternType, typename LocalLinkType> __global__ void kernelCalculatePolyakovLoops( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn, int timeDir = 0 )
	{
		typedef typename PatternType::PARAMTYPE::REALTYPE T;
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site(dim, nn );
		site.setLatticeIndex( index );

		WilsonLoop<PatternType,LocalLinkType> w( U, dim );

		data[index] = w.getReTracePolyakovLoopNormalized( site, timeDir );
	};
}

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> > > class WilsonLoopAverageCudaHelper
{
public:

	static void calculateWilsonLoops( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, int dir1, int length1, int dir2, int length2 )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		WilsonLoopAverageCudaHelperKernel::kernelCalculateWilsonLoops<PatternType,LocalLinkType><<<setup.getGridSize(),setup.getBlockSize()>>>( U, data, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ), dir1, length1, dir2, length2 );
		CUDA_LAST_ERROR( "kernelCalculateWilsonLoops" );
	}

	static void calculatePolyakovLoops( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, int timeDir = 0 )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		WilsonLoopAverageCudaHelperKernel::kernelCalculatePolyakovLoops<PatternType,LocalLinkType><<<setup.getGridSize(),setup.getBlockSize()>>>( U, data, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ), timeDir );
		CUDA_LAST_ERROR( "kernelCcalculatePolyakovLoopCorrelators" );
	}

	static void calculatePolyakovLoopCorrelators( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* data,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, int dirSeparation, int separation, int timeDir = 0 )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		WilsonLoopAverageCudaHelperKernel::kernelCalculatePolyakovLoopCorrelators<PatternType,LocalLinkType><<<setup.getGridSize(),setup.getBlockSize()>>>( U, data, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ), dirSeparation, separation, timeDir );
		CUDA_LAST_ERROR( "kernelCcalculatePolyakovLoopCorrelators" );
	}

};

}

#endif
