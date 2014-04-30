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

TEST_F( SU2RealFull, GetSU2Subgroup )
{
	store[1] = someValue;

	typename Real4<float>::VECTORTYPE subgroup = SUNRealFull<2,float>::getSU2Subgroup( store, 0, 1 );

	ASSERT_FLOAT_EQ( someValue, subgroup.y );
}

void fillSU2( float* store )
{
	store[0] = 0.1980;
	store[1] = 0.5189;
	store[2] = 0.4720;
	store[3] = 0.6847;
	store[4] = -0.4720;
	store[5] = 0.6847;
	store[6] = 0.1980;
	store[7] = -0.5189;
}

void assertRightSubgroup( float store[8] )
{
	float result[8] = {-0.8468152, 0.6493795,  0.5535316, 0.4073048,
			 -0.5535316, 0.4073048, -0.8468152, -0.6493795};
	for( int i = 0; i < 8; i++ )
	{
		ASSERT_FLOAT_EQ( result[i], store[i] );
	}
}

void assertLeftSubgroup( float store[8] )
{
	float result[8] = {-0.8468152, 0.3856970, 0.2254380, 0.8333088,
			 -0.2254380, 0.8333088, -0.8468152, -0.3856970};
	for( int i = 0; i < 8; i++ )
	{
		ASSERT_FLOAT_EQ( result[i], store[i] );
	}
}

TEST_F( SU2RealFull, RightSubgroupMult )
{
	fillSU2( store );

	typename Real4<float>::VECTORTYPE subgroup;
	subgroup.x = 0.7094;
	subgroup.y = 0.7547;
	subgroup.z = 0.2760;
	subgroup.w = 0.6797;

	SUNRealFull<2,float>::rightSubgroupMult( store, subgroup, 0, 1 );

	assertRightSubgroup( store );
}

TEST_F( SU2RealFull, LeftSubgroupMult )
{
	fillSU2( store );

	typename Real4<float>::VECTORTYPE subgroup;
	subgroup.x = 0.7094;
	subgroup.y = 0.7547;
	subgroup.z = 0.2760;
	subgroup.w = 0.6797;

	SUNRealFull<2,float>::leftSubgroupMult( store, subgroup, 0, 1 );

	assertLeftSubgroup( store );
}
