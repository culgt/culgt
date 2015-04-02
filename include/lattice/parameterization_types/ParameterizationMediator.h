/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_H_
#define PARAMETERIZATIONMEDIATOR_H_

#include "cudacommon/cuda_host_device.h"

namespace culgt
{

/*
 *
 */
template<typename ParamType1, typename ParamType2, typename LinkType1, typename LinkType2> class ParameterizationMediator
{
public:
	CUDA_HOST_DEVICE static void assign( LinkType1& l1, const LinkType2& l2 );
};

/**
 * Spezialization for same ParameterType but different LinkType
 */
template<typename LinkType1, typename LinkType2, typename ParamType> class ParameterizationMediator<ParamType,ParamType, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < ParamType::SIZE; i++ )
		{
			l1.set(i, l2.get(i));
		}
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_H_ */
