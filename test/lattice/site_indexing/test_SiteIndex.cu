/**
 * test_SiteIndex.cc
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "lattice/site_indexing/SiteNeighbourTableManager.h"
#include "lattice/LatticeDimension.h"
#include "testhelper_Site.h"
#include "cudatest/cuda_gtest_plugin.h"

using namespace ::testing;
using namespace culgt;


CUDA_TEST( ASiteIndex, GetCoord )
{
	LatticeDimension<4> dim(4,4,4,4);
	SiteIndex<4,NO_SPLIT> site( dim, DO_NOT_USE_NEIGHBOURS );
	site.setIndex( 123 );

	ASSERT_EQ( 1, site.getCoord(0) );
}
