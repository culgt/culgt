/**
 * test_SiteCoord.cc
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "lattice/site_indexing/SiteCoord.h"
#include "testhelper_Site.h"

using namespace ::testing;
using namespace culgt;


TEST(ASiteCoord, ShowSetGetUsage )
{
	// for SiteIndex get()/set() is trivial (see SiteIndex.hxx)
	// here it is non-trivial because we calculate coordinates from the index and
	// then again the index.
	int size[4] = {4,4,4,4};
	SiteCoord<4,NO_SPLIT> site( size );
	site.setIndex( 123 );

	ASSERT_EQ( 123, site.getIndex() );
}




class SiteCoordCompatibilityFullSplitNoSplit: public Test
{
public:
	static const int latticeIndex = 123;
	static const int size[4]; // size is 4^3 = 265
	SiteCoord<4,NO_SPLIT> siteNoSplit;
	SiteCoord<4,FULL_SPLIT> siteFullSplit;

	SiteCoordCompatibilityFullSplitNoSplit() : siteNoSplit(size), siteFullSplit(size)
	{
	}
};

const int SiteCoordCompatibilityFullSplitNoSplit::size[] = {4,4,4,4};

TEST_F( SiteCoordCompatibilityFullSplitNoSplit, CoordinatesAreEqual )
{
	siteNoSplit.setIndex( latticeIndex );
	siteFullSplit.setIndexFromNonParitySplitOrder( latticeIndex );

	ASSERT_THAT( true, coordinatesAreEqual( siteNoSplit, siteFullSplit ));
}




class SiteCoordCompatibilityTimesliceSplitNoSplit: public Test
{
public:
	static const int latticeIndex = 123;
	static const int size[4]; // size is 4^4 = 265
	SiteCoord<4,NO_SPLIT> siteNoSplit;
	SiteCoord<4,TIMESLICE_SPLIT> siteTimesliceSplit;

	SiteCoordCompatibilityTimesliceSplitNoSplit() : siteNoSplit(size), siteTimesliceSplit(size)
	{
	}
};

const int SiteCoordCompatibilityTimesliceSplitNoSplit::size[] = {4,4,4,4};

TEST_F( SiteCoordCompatibilityTimesliceSplitNoSplit, CoordinatesAreEqual )
{
	siteNoSplit.setIndex( latticeIndex );
	siteTimesliceSplit.setIndexFromNonParitySplitOrder( latticeIndex );

	ASSERT_THAT( true, coordinatesAreEqual( siteNoSplit, siteTimesliceSplit ));
}
