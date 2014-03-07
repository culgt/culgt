/**
 * test_LatticeDimension.cu
 *
 *  Created on: Mar 6, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "../cudatest/cuda_gtest_plugin.h"
#include "LatticeDimension.h"

using namespace culgt;
using namespace ::testing;


CUDA_TEST( ALatticeDimension, ReturnsSizeIn0DirectionFor1dimLattice )
{
	LatticeDimension<1> dim( 4 );

	ASSERT_EQ( 4, dim.getDimension(0) );
}

CUDA_TEST( ALatticeDimension, ReturnsSizeIn2DirectionFor4dimLattice )
{
	LatticeDimension<4> dim( 4,5,6,7 );

	ASSERT_EQ( 6, dim.getDimension(2) );
}

CUDA_TEST( ALatticeDimension, ReturnsLatticeSize )
{
	LatticeDimension<4> dim( 4,5,6,7 );

	ASSERT_EQ( 4*5*6*7, dim.getSize() );
}

CUDA_TEST( ALatticeDimension, SetSizeViaArray )
{
	int size[4] = {8,9,10,11};
	LatticeDimension<4> dim( size );

	ASSERT_EQ( 10, dim.getDimension(2) );
}

