/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "SUNRealFull.h"

using namespace culgt;
using namespace ::testing;

/**
 * SU3 is implicitly tested in test_LocalLink.cc
 */


//TEST( SU3RealFull, Identity )
//{
//	float store[18];
//
//	SUNRealFull<3,float>::identity( store );
//
//	ASSERT_FLOAT_EQ( 1.0, store[16] );
//}
//
//TEST( SU3RealFull, Det )
//{
//	float store[18];
//	SUNRealFull<3,float>::identity( store );
//
//	float result = SUNRealFull<3,float>::reDet( store );
//
//	ASSERT_FLOAT_EQ( 1.0, result );
//}
