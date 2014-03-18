/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "SU2Vector4.h"

using namespace culgt;
using namespace ::testing;


class ASU2Vector4: public Test
{
public:
	float4 store[1];
	float4 store2[1];
	static const float someValue = 1.42;
};

TEST_F( ASU2Vector4, Identity )
{
	SU2Vector4<float>::identity( store );

	ASSERT_FLOAT_EQ( 1.0, store[0].x );
}

TEST_F( ASU2Vector4, ReTrace )
{
	SU2Vector4<float>::identity( store );

	float result = SU2Vector4<float>::reTrace( store );

	ASSERT_FLOAT_EQ( 2.0, result );
}

TEST_F( ASU2Vector4, ReDet )
{
	SU2Vector4<float>::identity( store );

	float result = SU2Vector4<float>::reDet( store );

	ASSERT_FLOAT_EQ( 1.0, result );
}

TEST_F( ASU2Vector4, MultAssign )
{
	SU2Vector4<float>::identity( store );
	SU2Vector4<float>::identity( store2 );
	store2[0].x = someValue;

	SU2Vector4<float>::multAssign( store, store2 );

	ASSERT_FLOAT_EQ( someValue, store[0].x );
}

TEST_F( ASU2Vector4, Hermitian )
{
	SU2Vector4<float>::identity( store );
	store[0].w = 4.1;

	SU2Vector4<float>::hermitian( store );

	ASSERT_FLOAT_EQ( -4.1, store[0].w );
}
