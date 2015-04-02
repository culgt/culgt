/**
 * SiteNeighbourTableManager.h
 *
 *  Created on: Mar 17, 2014
 *      Author: vogt
 */

#ifndef SITENEIGHBOURTABLEMANAGER_H_
#define SITENEIGHBOURTABLEMANAGER_H_

#include <cstddef>
#include <map>
#include "lattice/LatticeDimension.h"
#include <iostream>

#ifdef __CUDACC__
#include "cudacommon/cuda_error.h"
#include "cudacommon/DeviceMemoryManager.h"
#endif


namespace culgt
{

template<typename SiteType> class SiteNeighbourTableManager
{
public:

	static lat_index_t* getHostPointer( LatticeDimension<SiteType::NDIM> dim )
	{
		if( !isAvailable( dim ) )
			generate(dim);

		return storeHost[dim];
	}

	static lat_index_t* getDevicePointer( LatticeDimension<SiteType::NDIM> dim )
	{
		if( !isAvailable( dim ) )
			generate(dim);

		return storeDevice[dim];
	}

	static bool isAvailable( const LatticeDimension<SiteType::NDIM> dim )
	{
		return (storeHost.count( dim )==1);
	}

	static void generate( const LatticeDimension<SiteType::NDIM> dim )
	{
		lat_index_t* nn = (lat_index_t*) malloc( dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ) );

		SiteType site( dim, nn );
		site.calculateNeighbourTable( nn );
		storeHost[dim] = nn;
#ifdef __CUDACC__
		lat_index_t* nnd;
		DeviceMemoryManager::malloc( &nnd, dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ), "Neighbour table" );
//		CUDA_SAFE_CALL( cudaMalloc( &nnd, dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ) ), "malloc neighbour table" );
		CUDA_SAFE_CALL( cudaMemcpy( nnd, nn, dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ), cudaMemcpyHostToDevice ), "memcpy neighbour table" );
		storeDevice[dim] = nnd;
#endif
	}

private:
	typedef std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> mymap;
	static mymap storeHost;
	static mymap storeDevice;
//	static SiteNeighbourTableManager<SiteType> manager;
};

template<typename SiteType> std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> SiteNeighbourTableManager<SiteType>::storeHost;
template<typename SiteType> std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> SiteNeighbourTableManager<SiteType>::storeDevice;

}

#endif /* SITENEIGHBOURTABLEMANAGER_H_ */


