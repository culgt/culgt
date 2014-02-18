/**
 * ParameterizationMediatorSU3_Real12_Real18_test.cc
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "ParameterizationMediatorSU3_Real12_Real18.h"

using namespace culgt;
using namespace ::testing;

class GetSetMock
{
public:
	float in[18];
	float out[12];
	float get( int i ) const
	{
		return in[i];
	}
	void set( int i, float val )
	{
		out[i] = val;
	}
};

class AParameterizationMediator_Real12_Real18: public Test
{
public:
	GetSetMock getset;
};

TEST_F( AParameterizationMediator_Real12_Real18, SpecializationCanBeCalled )
{
	ParameterizationMediator<SU3Real12<float>,SU3Real18<float>,GetSetMock,GetSetMock >::assign( getset, getset );
}

TEST_F( AParameterizationMediator_Real12_Real18, AssignCopiesData )
{
	getset.in[5] = 1.3;
	ParameterizationMediator<SU3Real12<float>,SU3Real18<float>,GetSetMock,GetSetMock >::assign( getset, getset );
	ASSERT_FLOAT_EQ( 1.3, getset.out[5] );
}
