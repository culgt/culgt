/**
 * test_SiteIndex.cc
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "SiteIndex.hxx"
#include "SiteCoord.hxx"
#include "../lattice/LatticeDimension.h"
#include <iostream>

using namespace std;
using namespace ::testing;
using namespace culgt;

TEST(ASiteIndex, ConvertSplitIndexToNoSplitIndex )
{
	LatticeDimension<4> dim( 4,6,7,3 );

	SiteIndex<4,FULL_SPLIT> siteIndexFullSplit( dim, DO_NOT_USE_NEIGHBOURS );

	SiteCoord<4,FULL_SPLIT> siteCoordFullSplit( dim );
	SiteCoord<4,NO_SPLIT> siteCoordNoSplit( dim );

	for( int t = 0; t < dim.getDimension( 0 ); t++ )
	for( int x = 0; x < dim.getDimension( 1 ); x++ )
	for( int y = 0; y < dim.getDimension( 2 ); y++ )
	for( int z = 0; z < dim.getDimension( 3 ); z++ )
	{
		siteCoordFullSplit[0] = t;
		siteCoordFullSplit[1] = x;
		siteCoordFullSplit[2] = y;
		siteCoordFullSplit[3] = z;

		siteCoordNoSplit[0] = t;
		siteCoordNoSplit[1] = x;
		siteCoordNoSplit[2] = y;
		siteCoordNoSplit[3] = z;

		siteIndexFullSplit.setIndex( siteCoordFullSplit.getIndex() );

//		cout << siteIndexFullSplit.getIndex() << endl;
		ASSERT_EQ( siteCoordNoSplit.getIndex(), siteIndexFullSplit.getIndexNonSplit() );
	}
}
