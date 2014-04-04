/**
 *  Created on: Apr 1, 2014
 *      Author: vogt
 */

#ifndef GPUPATTERNTIMESLICEPARITYPRIORITY_H_
#define GPUPATTERNTIMESLICEPARITYPRIORITY_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"

#include "StandardPattern.h"
#include "../LatticeDimension.h"
#include "../../cuLGT1legacy/SiteCoord.hxx"
#include <boost/mpl/assert.hpp>
#include "GPUPatternParityPriority.h"
#include "../../cuLGT1legacy/SiteIndex.hxx"

namespace culgt
{

template<typename Site, typename ParamType> class GPUPatternTimesliceParityPriority
{
private:
	BOOST_MPL_ASSERT_RELATION( Site::PARITYTYPE, ==, TIMESLICE_SPLIT );
public:
	typedef Site SITETYPE;
	typedef ParamType PARAMTYPE;
	typedef GPUPatternParityPriority<SiteIndex<Site::Ndim,FULL_SPLIT>, ParamType > TIMESLICE_PATTERNTYPE; // this should not be hard coded SiteIndex<Ndim,FULL_SPLIT> but any Site<Ndim,FULL_SPLIT> (depending on SITETYPE)

	CUDA_HOST_DEVICE static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		lat_index_t halfTimesliceSize = site.getSizeTimeslice()/2;
		lat_bool_t parity = site.getIndexTimeslice() / halfTimesliceSize;
		lat_index_t indexInParityInTimeslice = site.getIndexTimeslice() % halfTimesliceSize;

		return (((site[0]*2+parity)*SITETYPE::Ndim+mu)*ParamType::SIZE+paramIndex)*halfTimesliceSize+indexInParityInTimeslice;
	}

	CUDA_HOST_DEVICE static lat_array_index_t convertToStandardIndex( lat_array_index_t index, LatticeDimension<Site::Ndim> dim )
	{
		lat_index_t halfTimesliceSize = dim.getSizeTimeslice()/2;

		lat_index_t siteIndexInTimesliceInParity = index % halfTimesliceSize;

		index /= halfTimesliceSize;
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		index /= ParamType::SIZE;
		lat_dim_t mu = index % SITETYPE::Ndim;

		index /= SITETYPE::Ndim;
		lat_coord_t parity = index % 2;

		lat_coord_t t = index / 2;

		SiteCoord<Site::Ndim, NO_SPLIT> site( dim );
		site.setLatticeIndex( t*dim.getSizeTimeslice()+parity*halfTimesliceSize+siteIndexInTimesliceInParity);
		return StandardPattern<Site, ParamType>::getStandardIndex( site.getLatticeIndex(), mu, paramIndex );
	}
};

}

#endif /* GPUPATTERNTIMESLICEPARITYPRIORITY_H_ */
