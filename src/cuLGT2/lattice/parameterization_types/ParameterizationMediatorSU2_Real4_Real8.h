/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_REAL4_REAL8_H_
#define PARAMETERIZATIONMEDIATOR_REAL4_REAL8_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../ParameterizationMediator.h"
#include "SU2Real4.h"
#include "SUNRealFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SU2Real4<T>,SUNRealFull<2,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 4; i++ )
		{
			l1.set(i, l2.get(i));
		}
	}
};

template<typename LinkType1, typename LinkType2> class ParameterizationMediator<SU2Real4<float>,SUNRealFull<2,float>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 4; i++ )
		{
			l1.set(i, l2.get(i));
		}
	}
};

/**
 * Copy from 12 to 18 parameterization needs reconstruction of "third row"
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNRealFull<2,T>,SU2Real4<T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 4; i++ )
		{
			l1.set(i, l2.get(i));
		}

		l1.set(4, -l1.get(2) );
		l1.set(5,  l1.get(3) );
		l1.set(6,  l1.get(0) );
		l1.set(7, -l1.get(1) );
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_ */



