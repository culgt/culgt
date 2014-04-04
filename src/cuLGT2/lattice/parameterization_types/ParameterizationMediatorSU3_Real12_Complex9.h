/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_REAL12_COMPLEX9_H_
#define PARAMETERIZATIONMEDIATOR_REAL12_COMPLEX9_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../ParameterizationMediator.h"
#include "SU3Real12.h"
#include "SUNComplexFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SU3Real12<T>,SUNComplexFull<3,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 6; i++ )
		{
			typename SUNComplexFull<3,T>::TYPE temp;
			temp = l2.get( i );
			l1.set(i*2, temp.x);
			l1.set(i*2+1, temp.y);
		}
	}
};

//template<typename LinkType1, typename LinkType2> class ParameterizationMediator<SU3Real12<float>,SUNRealFull<3,float>, LinkType1, LinkType2>
//{
//public:
//	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
//	{
//		for( lat_group_index_t i = 0; i < 12; i++ )
//		{
//			l1.set(i, l2.get(i));
//		}
//	}
//};

/**
 * Copy from 12 to 18 parameterization needs reconstruction of "third row"
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNComplexFull<3,T>,SU3Real12<T>, LinkType1, LinkType2>
{
public:
	static CUDA_HOST_DEVICE inline lat_group_index_t getIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return i*3+j;
	}
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typedef typename SUNComplexFull<3,T>::TYPE COMPLEX;
		for( lat_group_index_t i = 0; i < 6; i++ )
		{
			COMPLEX temp;
			temp = COMPLEX( l2.get(2*i), l2.get(2*i+1) );
			l1.set( i, temp );
		}

		COMPLEX a = l1.get(getIndex(0,1)).conj() * l1.get(getIndex(1,2)).conj() - l1.get(getIndex(0,2)).conj() * l1.get(getIndex(1,1)).conj();
		COMPLEX b = l1.get(getIndex(0,2)).conj() * l1.get(getIndex(1,0)).conj() - l1.get(getIndex(0,0)).conj() * l1.get(getIndex(1,2)).conj();
		COMPLEX c = l1.get(getIndex(0,0)).conj() * l1.get(getIndex(1,1)).conj() - l1.get(getIndex(0,1)).conj() * l1.get(getIndex(1,0)).conj();

		T norm = a.abs_squared() + b.abs_squared() + c.abs_squared();
		norm = 1./sqrt( norm );

		l1.set( 6, a*norm );
		l1.set( 7, b*norm );
		l1.set( 8, c*norm );
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_ */



