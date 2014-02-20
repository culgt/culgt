/**
 * test_SiteIndex.cc
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "SiteIndex.hxx"
#include "testhelper_Site.h"

using namespace ::testing;

TEST(ASiteIndex, ShowSetGetUsage )
{
	int size[4] = {4,4,4,4};
	SiteIndex<4,NO_SPLIT> site( size );
	site.setLatticeIndex( 123 );

	ASSERT_EQ( 123, site.getLatticeIndex() );
}

class SiteIndexCompatibilityFullSplitNoSplit: public Test
{
public:
	static const int latticeIndex = 123;
	static const int size[4]; // size is 4^3 = 265
	SiteIndex<4,NO_SPLIT> siteNoSplit;
	SiteIndex<4,FULL_SPLIT> siteFullSplit;

	SiteIndexCompatibilityFullSplitNoSplit() : siteNoSplit(size), siteFullSplit(size)
	{
	}
};

const int  SiteIndexCompatibilityFullSplitNoSplit::size[] = {4,4,4,4};

// TODO: This is currently not implemented for SiteIndex
TEST_F( SiteIndexCompatibilityFullSplitNoSplit, DISABLED_IndexIsDifferent )
{
	siteNoSplit.setLatticeIndex( latticeIndex );
	siteFullSplit.setLatticeIndexFromNonParitySplitOrder( latticeIndex );

	ASSERT_THAT( true, coordinatesAreEqual( siteNoSplit, siteFullSplit ));
}
