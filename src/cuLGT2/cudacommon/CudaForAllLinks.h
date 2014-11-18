/**
 * ClassToExecute should be a class with operator()(LocalLink) defined.
 * This operator is called.
 *
 *  Created on: Apr 7, 2014
 *      Author: vogt
 */

#ifndef CUDAFORALLLINKS_H_
#define CUDAFORALLLINKS_H_
#include "../lattice/LatticeDimension.h"
#include "../lattice/KernelSetup.h"
#include "../cuLGT1legacy/SiteNeighbourTableManager.h"

namespace culgt
{


namespace CudaForAllLinksKernel
{
	template<typename PatternType, typename LocalLinkType, typename ClassToExecute> __global__ void execute( typename PatternType::PARAMTYPE::TYPE* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn )
	{
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site( dim, nn );
		site.setLatticeIndex( index );

		for( int mu = 0; mu < PatternType::SITETYPE::Ndim; mu++ )
		{
			GlobalLink<PatternType> glob( U, site, mu );
			LocalLinkType link;
			link = glob;
			ClassToExecute exec;
			exec( link );
			glob = link;
		}
	}

	template<typename PatternType, typename LocalLinkType, typename ClassToExecute> __global__ void executeSpatial( typename PatternType::PARAMTYPE::TYPE* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn )
	{
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site( dim, nn );
		site.setLatticeIndex( index );

		for( int mu = 1; mu < PatternType::SITETYPE::Ndim; mu++ )
		{
			GlobalLink<PatternType> glob( U, site, mu );
			LocalLinkType link;
			link = glob;
			ClassToExecute exec;
			exec( link );
			glob = link;
		}
	}

	template<typename PatternType, typename LocalLinkType, typename ClassToExecute> __global__ void execute( typename PatternType::PARAMTYPE::TYPE* U1, typename PatternType::PARAMTYPE::TYPE* U2, LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn )
	{
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site( dim, nn );
		site.setLatticeIndex( index );

		for( int mu = 0; mu < PatternType::SITETYPE::Ndim; mu++ )
		{
			GlobalLink<PatternType> glob1( U1, site, mu );
			LocalLinkType link1;
			link1 = glob1;
			GlobalLink<PatternType> glob2( U2, site, mu );
			LocalLinkType link2;
			link2 = glob2;
			ClassToExecute exec;
			exec( link1, link2 );
			glob1 = link1;
			glob2 = link2;
		}
	}
	template<typename PatternType, typename LocalLinkType, typename ClassToExecute> __global__ void executeSpatial( typename PatternType::PARAMTYPE::TYPE* U1, typename PatternType::PARAMTYPE::TYPE* U2, LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn )
	{
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site( dim, nn );
		site.setLatticeIndex( index );

		for( int mu = 1; mu < PatternType::SITETYPE::Ndim; mu++ )
		{
			GlobalLink<PatternType> glob1( U1, site, mu );
			LocalLinkType link1;
			link1 = glob1;
			GlobalLink<PatternType> glob2( U2, site, mu );
			LocalLinkType link2;
			link2 = glob2;
			ClassToExecute exec;
			exec( link1, link2 );
			glob1 = link1;
			glob2 = link2;
		}
	}
}


template<typename PatternType, typename LocalLinkType, typename ClassToExecute> class CudaForAllLinks
{
public:
	static void execute( typename PatternType::PARAMTYPE::TYPE* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		CudaForAllLinksKernel::execute<PatternType,LocalLinkType,ClassToExecute><<<setup.getGridSize(), setup.getBlockSize()>>>( U, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
	}

	static void executeSpatial( typename PatternType::PARAMTYPE::TYPE* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		CudaForAllLinksKernel::executeSpatial<PatternType,LocalLinkType,ClassToExecute><<<setup.getGridSize(), setup.getBlockSize()>>>( U, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
	}

	static void execute( typename PatternType::PARAMTYPE::TYPE* U1, typename PatternType::PARAMTYPE::TYPE* U2, LatticeDimension<PatternType::SITETYPE::Ndim> dim )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		CudaForAllLinksKernel::execute<PatternType,LocalLinkType,ClassToExecute><<<setup.getGridSize(), setup.getBlockSize()>>>( U1, U2, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
	}

	static void executeSpatial( typename PatternType::PARAMTYPE::TYPE* U1, typename PatternType::PARAMTYPE::TYPE* U2, LatticeDimension<PatternType::SITETYPE::Ndim> dim )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		CudaForAllLinksKernel::executeSpatial<PatternType,LocalLinkType,ClassToExecute><<<setup.getGridSize(), setup.getBlockSize()>>>( U1, U2, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
	}
};

}

#endif /* CUDAFORALLLINKS_H_ */
