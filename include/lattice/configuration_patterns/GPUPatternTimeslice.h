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

#endif /* GPUCOULOMBPATTERN_HXX_ */
