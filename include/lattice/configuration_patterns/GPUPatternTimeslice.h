/**
 */

#ifndef GPUPATTERNTIMESLICE_HXX_
#define GPUPATTERNTIMESLICE_HXX_

#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"

#include "StandardPattern.h"
#include "lattice/LatticeDimension.h"
#include "lattice/site_indexing/SiteCoord.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "GPUPattern.h"

namespace culgt
{

template<typename Site, typename ParamType> class GPUPatternTimeslice
{
public:
	typedef Site SITETYPE;
	typedef ParamType PARAMTYPE;
	typedef GPUPattern<SiteIndex<Site::NDIM,NO_SPLIT>, ParamType > TIMESLICE_PATTERNTYPE; // this should not be hard coded SiteIndex<Ndim,NO_SPLIT> but any Site<Ndim,NO_SPLIT> (depending on SITETYPE)

	CUDA_HOST_DEVICE static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		return ((site.getCoord(0)*SITETYPE::NDIM+mu)*ParamType::SIZE+paramIndex)*site.getSizeTimeslice()+site.getIndexTimeslice();
	}

	CUDA_HOST_DEVICE static lat_array_index_t convertToStandardIndex( lat_array_index_t index, LatticeDimension<Site::NDIM> dim )
	{
		lat_index_t siteIndexTimeslice = index % dim.getSizeTimeslice();

		index /= dim.getSizeTimeslice();
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		index /= ParamType::SIZE;
		lat_dim_t mu = index % SITETYPE::NDIM;

		lat_coord_t t = index / SITETYPE::NDIM;

		SiteCoord<Site::NDIM, NO_SPLIT> site( dim );
		site.setIndex( t*dim.getSizeTimeslice()+siteIndexTimeslice);
		return StandardPattern<Site, ParamType>::getStandardIndex( site.getIndex(), mu, paramIndex );
	}
};

}

//template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_index_t GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getSiteIndex( Site s )
//{
//	lat_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
//	return s.getLatticeIndexTimeslice()+s[0]*timesliceSize;
//}
//
//template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getLinkIndex( Site s, lat_dim_t mu )
//{
//	// TODO schreibs besser
//	lat_array_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
//	lat_array_index_t muSize = T_Nc*T_Nc*2*s.getLatticeSizeTimeslice();
//	return s.getLatticeIndexTimeslice()+mu*muSize+s[0]*timesliceSize;
//}
//
//template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
//{
//	lat_array_index_t timesliceSize = s.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
//	lat_array_index_t latSize = s.getLatticeSizeTimeslice();
//	return s.getLatticeIndexTimeslice() + latSize*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) )+s[0]*timesliceSize;
//}
//
////template< int T_Ndim, int T_Nc> int GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getIndex( int linkIndex, int i, int j, int c )
////{
////	return linkIndex +
////}
//
///**
// *
// */
//template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getUniqueIndex( Site s, lat_dim_t mu, lat_group_dim_t i, lat_group_dim_t j, bool c )
//{
//	assert(false);
//	return 0;
//}
//
///**
// * calculate the pattern index from unique index.
// */
//template<class Site, lat_dim_t T_Ndim, lat_group_dim_t T_Nc> lat_array_index_t GpuPatternTimeslice<Site, T_Ndim, T_Nc>::getIndexByUnique( lat_array_index_t uniqueIndex, lat_coord_t size[T_Ndim] )
//{
//
////	uniqueIndex /= site.getLatticeSize();
////	int mu = uniqueIndex /=
//
//
//	bool c = uniqueIndex % 2;
//	uniqueIndex /= 2;
//	lat_group_dim_t j = uniqueIndex % T_Nc;
//	uniqueIndex /= T_Nc;
//	lat_group_dim_t i = uniqueIndex % T_Nc;
//	uniqueIndex /= T_Nc;
//	lat_dim_t mu = uniqueIndex % T_Ndim;
//	uniqueIndex /= T_Ndim;
//	lat_array_index_t latticeIndex = uniqueIndex;
//
////	printf( "c %d\n", c );
////	printf( "j %d\n", j );
////	printf( "i %d\n", i );
////	printf( "mu %d\n", mu );
////	printf( "latind %d\n", latticeIndex );
//
//
//	Site site( size ); // parity must not be split for unique index
//	site.setLatticeIndexFromNonParitySplitOrder( latticeIndex );
//
//	lat_array_index_t timesliceSize = site.getLatticeSizeTimeslice()*T_Ndim*T_Nc*T_Nc*2;
//	lat_array_index_t latSize = site.getLatticeSizeTimeslice();
//	return site.getLatticeIndexTimeslice() + latSize*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) )+site[0]*timesliceSize;
//}

#endif /* GPUCOULOMBPATTERN_HXX_ */
