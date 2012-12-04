/*
 * GpuTimeslicePattern.hxx
 *
 *  Created on: Apr 13, 2012
 *      Author: vogt
 */

#ifndef STANDARDPATTERN_HXX_
#define STANDARDPATTERN_HXX_

#include "../cuda/cuda_host_device.h"
#include "../datatype/lattice_typedefs.h"

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> class StandardPattern
{
public:
	static const lat_group_dim_t Nc = T_Nc;
	static const lat_index_t Ndim = T_Ndim;
	CUDA_HOST_DEVICE static inline lat_index_t getSiteIndex( Site s );
	CUDA_HOST_DEVICE static inline lat_array_index_t getLinkIndex( Site s, lat_dim_t mu );
//	CUDA_HOST_DEVICE inline int getIndex( int linkIndex, int i, int j, int c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getUniqueIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] );
	CUDA_HOST_DEVICE static inline lat_array_index_t getUniqueIndex( lat_array_index_t index, lat_coord_t size[T_Ndim] );
};

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_index_t StandardPattern<Site, T_Ndim, T_Nc>::getSiteIndex( Site s )
{
	return  s.getLatticeIndex()*T_Ndim*T_Nc*T_Nc*2;
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t StandardPattern<Site, T_Ndim, T_Nc>::getLinkIndex( Site s, lat_dim_t mu )
{
	return (mu+T_Ndim*s.getLatticeIndex())*T_Nc*T_Nc*2;
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t StandardPattern<Site, T_Ndim, T_Nc>::getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
{
	return c+2*(j+T_Nc*(i+T_Nc*(mu+T_Ndim*s.getLatticeIndex())));
}

/**
 * Unique index is the same as the standard pattern index;
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t StandardPattern<Site, T_Ndim, T_Nc>::getUniqueIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
{
	return c+2*(j+T_Nc*(i+T_Nc*(mu+T_Ndim*s.getLatticeIndex())));
}

/**
 * Unique index is the same as the standard pattern index;
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t StandardPattern<Site, T_Ndim, T_Nc>::getUniqueIndex( lat_array_index_t index, lat_coord_t size[T_Ndim] )
{
	return index;
}

/**
 * Unique index is the same as the standard pattern index;
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t StandardPattern<Site, T_Ndim, T_Nc>::getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] )
{
	return uniqueIndex;
}

#endif /* STANDARDPATTERN_HXX_ */
