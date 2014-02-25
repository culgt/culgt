/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef LINK_H_
#define LINK_H_

#include "../cudacommon/cuda_host_device.h"

namespace culgt
{

/*
 *
 */
template<typename ParamType> class Link
{
public:
	CUDA_HOST_DEVICE virtual typename ParamType::TYPE get( lat_group_index_t i ) const = 0;
	CUDA_HOST_DEVICE virtual void set( lat_group_index_t i, typename ParamType::TYPE val ) = 0;
	CUDA_HOST_DEVICE virtual ~Link(){};
};

} /* namespace culgt */
#endif /* LINK_H_ */
