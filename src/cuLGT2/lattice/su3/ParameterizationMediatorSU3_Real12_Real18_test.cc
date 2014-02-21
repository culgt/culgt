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
	float data[18];
	float get( int i ) const
	{
		return data[i];
	}
	void set( int i, float val )
	{
		data[i] = val;
	}
};

class AParameterizationMediator_Real12_Real18: public Test
{
public:
	GetSetMock getset1;
	GetSetMock getset2;
};

TEST_F( AParameterizationMediator_Real12_Real18, SpecializationCanBeCalled )
{
	ParameterizationMediator<SU3Real12<float>,SU3Real18<float>,GetSetMock,GetSetMock >::assign( getset1, getset2 );
}

TEST_F( AParameterizationMediator_Real12_Real18, AssignCopiesData )
{
	getset2.data[5] = 1.3;
	ParameterizationMediator<SU3Real12<float>,SU3Real18<float>,GetSetMock,GetSetMock >::assign( getset1, getset2 );
	ASSERT_FLOAT_EQ( 1.3, getset1.data[5] );
}
