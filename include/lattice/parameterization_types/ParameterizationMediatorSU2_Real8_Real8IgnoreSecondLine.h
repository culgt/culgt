/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_REAL8_REAL8_IGNORESECONDLINE_H_
#define PARAMETERIZATIONMEDIATOR_REAL8_REAL8_IGNORESECONDLINE_H_

#include "../../cudacommon/cuda_host_device.h"
#include "../../common/culgt_typedefs.h"
#include "../ParameterizationMediator.h"
#include "SUNRealFull.h"
#include "SU2RealFullIgnoreSecondLine.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SUNRealFull<2,T>,SU2RealFullIgnoreSecondLine<T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		l1.set( 0, l2.get(0) );
		l1.set( 1, l2.get(1) );
		l1.set( 2, l2.get(2) );
		l1.set( 3, l2.get(3) );

		l1.set( 4, -l1.get(2) );
		l1.set( 5,  l1.get(3) );
		l1.set( 6,  l1.get(0) );
		l1.set( 7, -l1.get(1) );
	}
};

/**
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SU2RealFullIgnoreSecondLine<T>,SUNRealFull<2,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		l1.set( 0, l2.get( 0 ) );
		l1.set( 1, l2.get( 1 ) );
		l1.set( 2, l2.get( 2 ) );
		l1.set( 3, l2.get( 3 ) );
	}
};

} /* namespace culgt */

#endif



