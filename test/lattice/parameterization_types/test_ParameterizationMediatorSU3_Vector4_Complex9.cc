/**
 * ParameterizationMediatorSU3_Real12_Real18_test.cc
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Complex9.h"

using namespace culgt;
using namespace ::testing;

class GetSetMockComplex
{
public:
	Complex<float> data[9];
	Complex<float> get( int i ) const
	{
		return data[i];
	}
	void set( int i, Complex<float> val )
	{
		data[i] = val;
	}
};
class GetSetMockVector4
{
public:
	float4 data[3];
	float4 get( int i ) const
	{
		return data[i];
	}
	void set( int i, float4 val )
	{
		data[i] = val;
	}
};

class AParameterizationMediator_Vector4_Complex9: public Test
{
public:
	GetSetMockComplex getset1;
	GetSetMockVector4 getset2;
};

TEST_F( AParameterizationMediator_Vector4_Complex9, SpecializationCanBeCalled )
{
	ParameterizationMediator<SUNComplexFull<3,float>,SU3Vector4<float>,GetSetMockComplex,GetSetMockVector4 >::assign( getset1, getset2 );
}

TEST_F( AParameterizationMediator_Vector4_Complex9, AssignCopiesData )
{
	float4 a;
	a.x = 1.;
	a.y = 0.;
	a.z = 0.;
	a.w = 0.;
	float4 b;
	b.x = 0.;
	b.y = 0.;
	b.z = 0.;
	b.w = 0.;
	float4 c;
	c.x = 1.;
	c.y = 0.;
	c.z = 0.;
	c.w = 0.;
	getset2.data[0] = a;
	getset2.data[1] = b;
	getset2.data[2] = c;

	ParameterizationMediator<SUNComplexFull<3,float>,SU3Vector4<float>,GetSetMockComplex,GetSetMockVector4 >::assign( getset1, getset2 );

	ASSERT_FLOAT_EQ( 1., getset1.data[9].x );
}
