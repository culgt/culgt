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


class SU2RealFull: public Test
{
public:
	float store[8];
	float store2[8];
	const float someValue = 1.42;
};

TEST_F( SU2RealFull, Identity )
{
	SUNRealFull<2,float>::identity( store );

	ASSERT_FLOAT_EQ( 1.0, store[6] );
}

TEST_F( SU2RealFull, ReTrace )
{
	SUNRealFull<2,float>::identity( store );

	float result = SUNRealFull<2,float>::reTrace( store );

	ASSERT_FLOAT_EQ( 2.0, result );
}

TEST_F( SU2RealFull, ReDet )
{
	SUNRealFull<2,float>::identity( store );

	float result = SUNRealFull<2,float>::reDet( store );

	ASSERT_FLOAT_EQ( 1.0, result );
}

TEST_F( SU2RealFull, MultAssign )
{
	SUNRealFull<2,float>::identity( store );
	SUNRealFull<2,float>::identity( store2 );
	store2[0] = someValue;

	SUNRealFull<2,float>::multAssign( store, store2 );

	ASSERT_FLOAT_EQ( someValue, store[0] );
}

TEST_F( SU2RealFull, Hermitian )
{
	SUNRealFull<2,float>::identity( store );
	store[3] = someValue;

	SUNRealFull<2,float>::hermitian( store );

	ASSERT_FLOAT_EQ( -someValue, store[5] );
}
