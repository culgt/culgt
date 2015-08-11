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
#include <iostream>
#include "lattice/LatticeDimension.h"

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
		if( !isAvailableOnHost( dim ) )
			generateOnHost(dim);

		return storeHost[dim];
	}

#ifdef __CUDACC__
	static lat_index_t* getDevicePointer( LatticeDimension<SiteType::NDIM> dim )
	{
		if( !isAvailableOnDevice( dim ) )
		{
			if( !isAvailableOnHost( dim ) )
			{
				generateOnHost( dim );
			}
			copyToDevice( dim );
		}

		return storeDevice[dim];
	}
#endif

	static bool isAvailableOnHost( const LatticeDimension<SiteType::NDIM> dim )
	{
		return (storeHost.count( dim )==1);
	}

	static bool isAvailableOnDevice( const LatticeDimension<SiteType::NDIM> dim )
	{
		return (storeDevice.count( dim )==1);
	}

	static void generateOnHost( const LatticeDimension<SiteType::NDIM> dim )
	{
		lat_index_t* nn = new lat_index_t[dim.getSize()*SiteType::NDIM*2];

		SiteType site( dim, nn );
		site.calculateNeighbourTable( nn );
		storeHost[dim] = nn;
	}

#ifdef __CUDACC__
	static void copyToDevice( const LatticeDimension<SiteType::NDIM> dim )
	{
		lat_index_t* nnd;
		DeviceMemoryManager::malloc( &nnd, dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ), "Neighbour table" );
		CUDA_SAFE_CALL( cudaMemcpy( nnd, storeHost[dim], dim.getSize()*SiteType::NDIM*2*sizeof( lat_index_t ), cudaMemcpyHostToDevice ), "memcpy neighbour table" );
		storeDevice[dim] = nnd;
	}
#endif

private:
	typedef std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> mymap;
	static mymap storeHost;
	static mymap storeDevice;
};

template<typename SiteType> std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> SiteNeighbourTableManager<SiteType>::storeHost;
template<typename SiteType> std::map<LatticeDimension<SiteType::NDIM>,lat_index_t*,LatticeDimensionCompare> SiteNeighbourTableManager<SiteType>::storeDevice;

}

#endif /* SITENEIGHBOURTABLEMANAGER_H_ */


