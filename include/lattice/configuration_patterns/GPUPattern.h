/**
 * GPUPattern.h
 *
 * The pattern is as follows (slowest running index first):
 * mu, paramIndex, site (where paramIndex for float18 parameterization is a value from 0..17)
 *
 * TODO Test difference to paramIndex,mu,site
 *
 *  Created on: Feb 19, 2014
 *      Author: vogt
 */

#ifndef GPUPATTERN_H_
#define GPUPATTERN_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../LatticeDimension.h"

#include "StandardPattern.h"
namespace culgt
{


template<typename Site, typename ParamType> class GPUPattern
{
public:
	typedef Site SITETYPE;
	typedef ParamType PARAMTYPE;

	CUDA_HOST_DEVICE static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		return (mu*ParamType::SIZE+paramIndex)*site.getSize()+site.getIndex();
	}

	/**
	 * For inter-pattern compatibility every pattern should implement a function that converts
	 * its own index to the standard index (which is the index of StandardPattern).
	 * @param index
	 * @param latSize
	 * @return
	 */
	CUDA_HOST_DEVICE static lat_array_index_t convertToStandardIndex( lat_array_index_t index, LatticeDimension<Site::Ndim> dim )
	{
		lat_index_t siteIndex = index % dim.getSize();

		index /= dim.getSize();
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		lat_dim_t mu = index / ParamType::SIZE;

		return StandardPattern<Site, ParamType>::getStandardIndex( siteIndex, mu, paramIndex );
	}
};

}

#endif /* GPUPATTERN_H_ */
