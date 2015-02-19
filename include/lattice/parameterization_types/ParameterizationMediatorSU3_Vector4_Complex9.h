/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_SU3_VECTOR4_COMPLEX9_H_
#define PARAMETERIZATIONMEDIATOR_SU3_VECTOR4_COMPLEX9_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../ParameterizationMediator.h"
#include "SU3Vector4.h"
#include "SUNComplexFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SU3Vector4<T>,SUNComplexFull<3,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typename SU3Vector4<T>::TYPE temp;
		typename SUNComplexFull<3,T>::TYPE complextemp;

		complextemp = l2.get(0);
		temp.x = complextemp.x;
		temp.y = complextemp.y;
		complextemp = l2.get(1);
		temp.z = complextemp.x;
		temp.w = complextemp.y;
		l1.set( 0, temp );

		complextemp = l2.get(2);
		temp.x = complextemp.x;
		temp.y = complextemp.y;
		complextemp = l2.get(3);
		temp.z = complextemp.x;
		temp.w = complextemp.y;
		l1.set( 1, temp );

		complextemp = l2.get(4);
		temp.x = complextemp.x;
		temp.y = complextemp.y;
		complextemp = l2.get(5);
		temp.z = complextemp.x;
		temp.w = complextemp.y;
		l1.set( 2, temp );
	}
};

//template<typename LinkType1, typename LinkType2> class ParameterizationMediator<SU2Real4<float>,SUNRealFull<2,float>, LinkType1, LinkType2>
//{
//public:
//	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
//	{
//		for( lat_group_index_t i = 0; i < 4; i++ )
//		{
//			l1.set(i, l2.get(i));
//		}
//	}
//};

/**
 * Copy from 12 to 18 parameterization needs reconstruction of "third row"
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNComplexFull<3,T>,SU3Vector4<T>, LinkType1, LinkType2>
{
public:
	static CUDA_HOST_DEVICE inline lat_group_index_t getIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return i*3+j;
	}
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typedef typename SUNComplexFull<3,T>::TYPE COMPLEX;

		typename SU2Vector4<T>::TYPE temp = l2.get(0);
		l1.set( 0, COMPLEX( temp.x, temp.y ) );
		l1.set( 1, COMPLEX( temp.z, temp.w ) );
		temp = l2.get(1);
		l1.set( 2, COMPLEX( temp.x, temp.y ) );
		l1.set( 3, COMPLEX( temp.z, temp.w ) );
		temp = l2.get(2);
		l1.set( 4, COMPLEX( temp.x, temp.y ) );
		l1.set( 5, COMPLEX( temp.z, temp.w ) );

		COMPLEX a = l1.get(getIndex(0,1)).conj() * l1.get(getIndex(1,2)).conj() - l1.get(getIndex(0,2)).conj() * l1.get(getIndex(1,1)).conj();
		COMPLEX b = l1.get(getIndex(0,2)).conj() * l1.get(getIndex(1,0)).conj() - l1.get(getIndex(0,0)).conj() * l1.get(getIndex(1,2)).conj();
		COMPLEX c = l1.get(getIndex(0,0)).conj() * l1.get(getIndex(1,1)).conj() - l1.get(getIndex(0,1)).conj() * l1.get(getIndex(1,0)).conj();

		T norm = a.abs_squared() + b.abs_squared() + c.abs_squared();
		norm = 1./::sqrt( norm );

		l1.set( 6, a*norm );
		l1.set( 7, b*norm );
		l1.set( 8, c*norm );
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_ */



