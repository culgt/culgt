/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU2_Real4_Real8.h"

using namespace culgt;
using namespace ::testing;

class GetSetMock
{
public:
	float data[8];
	float get( int i ) const
	{
		return data[i];
	}
	void set( int i, float val )
	{
		data[i] = val;
	}
};

class AParameterizationMediator_Real4_Real8: public Test
{
public:
	GetSetMock getset1;
	GetSetMock getset2;
};

TEST_F( AParameterizationMediator_Real4_Real8, AssignCopiesData )
{
	getset2.data[1] = 1.3;
	ParameterizationMediator<SUNRealFull<2,float>,SU2Real4<float>,GetSetMock,GetSetMock >::assign( getset1, getset2 );
	ASSERT_FLOAT_EQ( -1.3, getset1.data[7] );
}
