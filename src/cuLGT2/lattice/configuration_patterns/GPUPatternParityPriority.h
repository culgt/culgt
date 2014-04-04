/**
 * GPUPatternParityPriority.h
 *
 * The pattern is as follows (slowest running index first):
 * parity, mu, paramIndex, siteInParity (where paramIndex for float18 parameterization is a value from 0..17)
 *
 * IMPORTANT: Use a Site with FULL_SPLIT parity type.
 */

#ifndef GPUPATTERNPARITYPRIORITY_H_
#define GPUPATTERNPARITYPRIORITY_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "StandardPattern.h"
#include "../LatticeDimension.h"
#include <boost/mpl/assert.hpp>


namespace culgt
{


template<typename Site, typename ParamType> class GPUPatternParityPriority
{
private:
	BOOST_MPL_ASSERT_RELATION( Site::PARITYTYPE, ==, FULL_SPLIT ); // verify that a FULL_SPLIT site is used in combination with this pattern
public:
	typedef Site SITETYPE;
	typedef ParamType PARAMTYPE;

	CUDA_HOST_DEVICE static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		lat_index_t halfLatSize = site.getSize()/2;
		lat_bool_t parity = site.getIndex() / halfLatSize;
		lat_index_t indexInParity = site.getIndex() % halfLatSize;
		return ((parity*Site::Ndim+mu)*ParamType::SIZE+paramIndex)*halfLatSize+indexInParity;
	}

	/**
	 * For intern-pattern compatibility every pattern should implement a function that converts
	 * its own index to the standard index (which is the index of StandardPattern).
	 * @param index
	 * @param latSize
	 * @return
	 */
	CUDA_HOST_DEVICE static lat_array_index_t convertToStandardIndex( lat_array_index_t index, LatticeDimension<Site::Ndim> dim )
	{
		lat_index_t halfLatSize = dim.getSize()/2;

		lat_index_t siteIndexInParity = index % halfLatSize;

		index /= halfLatSize;
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		index /= ParamType::SIZE;
		lat_dim_t mu = index % Site::Ndim;

		lat_bool_t parity = index / Site::Ndim;

		lat_index_t siteIndex = siteIndexInParity + parity*halfLatSize;

		return StandardPattern<Site, ParamType>::getStandardIndex( siteIndex, mu, paramIndex );
	}
};

}

#endif /* GPUPATTERNPARITYPRIORITY_H_ */

