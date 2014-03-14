/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "SU3Real12.h"

using namespace culgt;
using namespace ::testing;

TEST( ASU3Real12, Identity )
{
	float store[12];

	SU3Real12<float>::identity( store );

	ASSERT_FLOAT_EQ( 1.0, store[8] );
	ASSERT_FLOAT_EQ( 1.0, store[0] );
}

