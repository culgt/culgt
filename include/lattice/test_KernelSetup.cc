/**
 * test_KernelSetup.cc
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "KernelSetup.h"
#include "LatticeDimension.h"

using namespace culgt;

TEST( AKernelSetup, SimpleSetup )
{
	LatticeDimension<4> dim( 4,4,4,4 );
	KernelSetup<4> kernel( dim );

	ASSERT_EQ( 8, kernel.getGridSize() );
}

TEST( AKernelSetup, SimpleSetupWithParity )
{
	LatticeDimension<4> dim( 4,4,4,4 );
	KernelSetup<4> kernel( dim, true );

	ASSERT_EQ( 4, kernel.getGridSize() );
}

TEST( AKernelSetup, NonDividableSetup )
{
	LatticeDimension<4> dim( 3,3,3,3 ); // 3^4 == 81, 81/32 = 2, 81%32=17
	KernelSetup<4> kernel( dim );

	ASSERT_EQ( 3, kernel.getGridSize() );
}
