/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_VECTOR4_COMPLEX4_H_
#define PARAMETERIZATIONMEDIATOR_VECTOR4_COMPLEX4_H_

#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"
#include "ParameterizationMediator.h"
#include "SU2Vector4.h"
#include "SUNComplexFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SU2Vector4<T>,SUNComplexFull<2,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typename SU2Vector4<T>::TYPE temp;
		Complex<T> templ2 = l2.get(0);
		temp.x = templ2.x;
		temp.y = templ2.y;
		templ2 = l2.get(1);
		temp.z = templ2.x;
		temp.w = templ2.y;
		l1.set( 0, temp );
	}
};


/**
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNComplexFull<2,T>,SU2Vector4<T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		typename SU2Vector4<T>::TYPE temp = l2.get(0);

		l1.set( 0, Complex<T>( temp.x, temp.y ) );
		l1.set( 1, Complex<T>( temp.z, temp.w ) );

		l1.set( 2, Complex<T>( -temp.z, temp.w ) );
		l1.set( 3, Complex<T>( temp.x, -temp.y ) );
	}
};

} /* namespace culgt */

#endif



