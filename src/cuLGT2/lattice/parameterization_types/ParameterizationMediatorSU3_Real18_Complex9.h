/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_REAL18_COMPLEX9_H_
#define PARAMETERIZATIONMEDIATOR_REAL18_COMPLEX9_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../ParameterizationMediator.h"
#include "SUNComplexFull.h"
#include "SUNRealFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SUNRealFull<3,T>,SUNComplexFull<3,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 9; i++ )
		{
			typename SUNComplexFull<3,T>::TYPE temp;
			temp = l2.get( i );
			l1.set(i*2, temp.x);
			l1.set(i*2+1, temp.y);
		}
	}
};

template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNComplexFull<3,T>,SUNRealFull<3,T>, LinkType1, LinkType2>
{
public:
	static CUDA_HOST_DEVICE inline lat_group_index_t getIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return i*3+j;
	}
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typedef typename SUNComplexFull<3,T>::TYPE COMPLEX;
		for( lat_group_index_t i = 0; i < 9; i++ )
		{
			COMPLEX temp;
			temp = COMPLEX( l2.get(2*i), l2.get(2*i+1) );
			l1.set( i, temp );
		}
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_ */



