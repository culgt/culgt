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

using namespace ::testing;
using culgt::SiteNeighbourTableManager;
using culgt::LatticeDimension;
using namespace culgt;

TEST(ASiteIndex, ShowSetGetUsage )
{
	int size[4] = {4,4,4,4};
	SiteIndex<4,NO_SPLIT> site( size, DO_NOT_USE_NEIGHBOURS );
	site.setIndex( 123 );

	ASSERT_EQ( 123, site.getIndex() );
}

TEST(ASiteIndex, InitializeWithLatticeDimensionAndNoNeighbours )
{
	LatticeDimension<4> dim(4,4,4,4);
	SiteIndex<4,NO_SPLIT> site( dim, DO_NOT_USE_NEIGHBOURS );
	site.setIndex( 123 );

	ASSERT_EQ( 123, site.getIndex() );
}

class SiteIndexCompatibilityFullSplitNoSplit: public Test
{
public:
	static const int latticeIndex = 123;
	static const int size[4]; // size is 4^3 = 265
	SiteIndex<4,NO_SPLIT> siteNoSplit;
	SiteIndex<4,FULL_SPLIT> siteFullSplit;

	SiteIndexCompatibilityFullSplitNoSplit() : siteNoSplit(size, DO_NOT_USE_NEIGHBOURS), siteFullSplit(size, DO_NOT_USE_NEIGHBOURS)
	{
	}
};

const int  SiteIndexCompatibilityFullSplitNoSplit::size[] = {4,4,4,4};

// TODO: To compare the coordinates we need FullSplit getCoord()
TEST_F( SiteIndexCompatibilityFullSplitNoSplit, DISABLED_IndexIsDifferent )
{
	siteNoSplit.setIndex( latticeIndex );
	siteFullSplit.setIndexFromNonParitySplitOrder( latticeIndex );

	ASSERT_THAT( true, coordinatesAreEqual( siteNoSplit, siteFullSplit ));
}

class SiteIndexCompatibilityFullSplitTimesliceSplit: public Test
{
public:
	const LatticeDimension<4> dim; // size is 4^3 = 265
	const LatticeDimension<4> dimTimeslice; // size is 4^3 = 265
	SiteIndex<4,TIMESLICE_SPLIT> siteTimesliceSplit;
	SiteIndex<4,FULL_SPLIT> siteFullSplitInTimeslice;

	SiteIndexCompatibilityFullSplitTimesliceSplit() : dim(4,4,4,4), dimTimeslice(1,4,4,4), siteTimesliceSplit(dim,SiteNeighbourTableManager<SiteIndex<4,TIMESLICE_SPLIT> >::getHostPointer( dim )), siteFullSplitInTimeslice(dimTimeslice,SiteNeighbourTableManager<SiteIndex<4,FULL_SPLIT> >::getHostPointer( dimTimeslice ))
	{
	}
};

TEST_F( SiteIndexCompatibilityFullSplitTimesliceSplit, CheckCompatibiltyInTimesliceSublattice )
{
	int site[4] = {2,3,3,3};
	int t = 2;
	int siteInTimeslice[4] = {0,3,3,3};
	siteTimesliceSplit.setIndex( siteTimesliceSplit.getIndex( site ) );
	siteFullSplitInTimeslice.setIndex( siteFullSplitInTimeslice.getIndex( siteInTimeslice ) );

	ASSERT_EQ(  siteTimesliceSplit.getIndex(), siteFullSplitInTimeslice.getIndex() + t*dimTimeslice.getSize() );
}

TEST_F( SiteIndexCompatibilityFullSplitTimesliceSplit, CheckNeighbourCompatibiltyInTimesliceSublattice )
{
	int site[4] = {2,3,3,3};
	int t = 2;
	int siteInTimeslice[4] = {0,3,3,3};
	siteTimesliceSplit.setIndex( siteTimesliceSplit.getIndex( site ) );
	siteFullSplitInTimeslice.setIndex( siteFullSplitInTimeslice.getIndex( siteInTimeslice ) );

	siteTimesliceSplit.setNeighbour( 1, true );
	siteFullSplitInTimeslice.setNeighbour( 1, true );

	ASSERT_EQ(  siteTimesliceSplit.getIndex(), siteFullSplitInTimeslice.getIndex() + t*dimTimeslice.getSize() );
}
