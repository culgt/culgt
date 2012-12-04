/*
 * GpuTimeslicePattern.hxx
 *
 *  Created on: Apr 13, 2012
 *      Author: vogt
 */

#ifndef GPUCOULOMBPATTERN_HXX_
#define GPUCOULOMBPATTERN_HXX_

#include <assert.h>
#include "../cuda/cuda_host_device.h"
#include "../datatype/lattice_typedefs.h"

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> class GpuCoulombPattern
{
public:
	static const lat_group_dim_t Nc = T_Nc; // TODO static is shit. do it via constructor and non static
	static const lat_index_t Ndim = T_Ndim;
	CUDA_HOST_DEVICE static inline lat_index_t getSiteIndex( Site s );
	CUDA_HOST_DEVICE static inline lat_array_index_t getLinkIndex( Site s, lat_dim_t mu );
//	CUDA_HOST_DEVICE static inline int getIndex( int linkIndex, int i, int j, int c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getUniqueIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c );
	CUDA_HOST_DEVICE static inline lat_array_index_t getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] );
};

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_index_t GpuCoulombPattern<Site, T_Ndim, T_Nc>::getSiteIndex( Site s )
{
	lat_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
	return s.getLatticeIndexTimeslice()+s[0]*timesliceSize;
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuCoulombPattern<Site, T_Ndim, T_Nc>::getLinkIndex( Site s, lat_dim_t mu )
{
	// TODO schreibs besser
	lat_array_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
	lat_array_index_t muSize = T_Nc*T_Nc*2*s.getLatticeSizeTimeslice();
	return s.getLatticeIndexTimeslice()+mu*muSize+s[0]*timesliceSize;
}

template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuCoulombPattern<Site, T_Ndim, T_Nc>::getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
{
	lat_array_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
	lat_array_index_t latSize = s.getLatticeSizeTimeslice();
	return s.getLatticeIndexTimeslice() + latSize*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) )+s[0]*timesliceSize;
}

//template< int T_Ndim, int T_Nc> int GpuCoulombPattern<Site, T_Ndim, T_Nc>::getIndex( int linkIndex, int i, int j, int c )
//{
//	return linkIndex +
//}

/**
 *
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuCoulombPattern<Site, T_Ndim, T_Nc>::getUniqueIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
{
	assert(false);
	return 0;
}

/**
 * calculate the pattern index from unique index.
 */
template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuCoulombPattern<Site, T_Ndim, T_Nc>::getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] )
{

//	uniqueIndex /= site.getLatticeSize();
//	int mu = uniqueIndex /=


	bool c = uniqueIndex % 2;
	uniqueIndex /= 2;
	lat_group_dim_t j = uniqueIndex % T_Nc;
	uniqueIndex /= T_Nc;
	lat_group_dim_t i = uniqueIndex % T_Nc;
	uniqueIndex /= T_Nc;
	lat_dim_t mu = uniqueIndex % T_Ndim;
	uniqueIndex /= T_Ndim;
	lat_array_index_t latticeIndex = uniqueIndex;

//	printf( "c %d\n", c );
//	printf( "j %d\n", j );
//	printf( "i %d\n", i );
//	printf( "mu %d\n", mu );
//	printf( "latind %d\n", latticeIndex );


	Site site( size ); // parity must not be split for unique index
	site.setLatticeIndexFromNonParitySplitOrder( latticeIndex );

	lat_array_index_t timesliceSize = site.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
	lat_array_index_t latSize = site.getLatticeSizeTimeslice();
	return site.getLatticeIndexTimeslice() + latSize*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) )+site[0]*timesliceSize;
}

#endif /* GPUCOULOMBPATTERN_HXX_ */
