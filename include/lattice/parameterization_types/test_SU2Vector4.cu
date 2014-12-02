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
	const static float someValue = 1.42;
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

TEST_F( ASU2Vector4, MultAssignScalar )
{
	SU2Vector4<float>::identity( store );
	const float someFactor = 0.5;
	store[0].z = someValue;

	SU2Vector4<float>::multAssignScalar( store, someFactor );

	ASSERT_FLOAT_EQ( someValue*someFactor, store[0].z );
}

TEST_F( ASU2Vector4, AddAssign )
{
	SU2Vector4<float>::identity( store );
	SU2Vector4<float>::identity( store2 );
	store2[0].x = someValue;

	SU2Vector4<float>::addAssign( store, store2 );

	ASSERT_FLOAT_EQ( someValue+1., store[0].x );
}

TEST_F( ASU2Vector4, SubtractAssign )
{
	SU2Vector4<float>::identity( store );
	SU2Vector4<float>::identity( store2 );
	store2[0].x = someValue;

	SU2Vector4<float>::subtractAssign( store, store2 );

	ASSERT_FLOAT_EQ( 1.-someValue, store[0].x );
}

TEST_F( ASU2Vector4, FrobeniusNorm )
{
	SU2Vector4<float>::zero( store );
	store[0].x = someValue;

	ASSERT_FLOAT_EQ( 2.*(someValue*someValue), SU2Vector4<float>::normFrobeniusSquared( store ) );
}

TEST_F( ASU2Vector4, Hermitian )
{
	SU2Vector4<float>::identity( store );
	store[0].w = 4.1;

	SU2Vector4<float>::hermitian( store );

	ASSERT_FLOAT_EQ( -4.1, store[0].w );
}

TEST_F( ASU2Vector4, Reproject )
{
	store[0].x = .25;
	store[0].y = .25;
	store[0].z = .25;
	store[0].w = .25;

	SU2Vector4<float>::reproject( store );

	ASSERT_FLOAT_EQ( .5, store[0].w );
}


class ASpecialSU2Vector4: public Test
{
public:
	float4 A[1];
	float4 B;

	float4 AB;
	float4 BA;

	const static float someValue = 1.42;

	ASpecialSU2Vector4()
	{
		A[0].x = 0.8003;
		A[0].y = 0.4218;
		A[0].z = 0.1419;
		A[0].w = 0.9157;

		B.x = 0.8491;
		B.y = 0.6787;
		B.z = 0.9340;
		B.w = 0.7577;

		AB.x = -0.4331014;
		AB.y = 0.1535678;
		AB.z = 1.1698552;
		AB.w = 1.6815619;

		BA.x = -0.4331014;
		BA.y =  1.6490601;
		BA.z = 0.5660797;
		BA.w = 1.0862546;
	}

protected:
	void assertFloatEq( float4 expect, float4 actual )
	{
		ASSERT_FLOAT_EQ( expect.x, actual.x );
		ASSERT_FLOAT_EQ( expect.y, actual.y );
		ASSERT_FLOAT_EQ( expect.z, actual.z );
		ASSERT_FLOAT_EQ( expect.w, actual.w );
	}

};

TEST_F( ASpecialSU2Vector4, RightSubgroupMult )
{
	SU2Vector4<float>::rightSubgroupMult( A, B, 0, 0 );

	assertFloatEq( AB, A[0] );
}

TEST_F( ASpecialSU2Vector4, LeftSubgroupMult )
{
	SU2Vector4<float>::leftSubgroupMult( A, B, 0, 0 );

	assertFloatEq( BA, A[0] );
}
