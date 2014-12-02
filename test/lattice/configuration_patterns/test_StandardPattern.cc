#include "gmock/gmock.h"
#include "lattice/configuration_patterns/StandardPattern.h"
#include "testhelper_PatternMocks.h"

using namespace culgt;
using namespace ::testing;

TEST( StandardPatternTest, GetIndexReturnsParamIndexForSiteZeroAndMuZero )
{
	SiteTypeMock<> mySite(0);
	const int myParamIndex = 4;

	int result = StandardPattern<SiteTypeMock<>, ParamTypeMock<1> >::getIndex( mySite, 0, myParamIndex );

	ASSERT_EQ( myParamIndex, result );
}

TEST( StandardPatternTest, GetIndexReturnsCorrectIndexForSiteZeroWithMockParamType )
{
	SiteTypeMock<> mySite( 0 );
	const int myParamSize = 18;
	const int myParamIndex = 4;
	const int myMu = 1;
	int expect = (myMu)*myParamSize+myParamIndex;

	int result = StandardPattern<SiteTypeMock<>, ParamTypeMock<myParamSize> >::getIndex( mySite, myMu, myParamIndex );

	ASSERT_EQ( expect, result );
}

TEST( StandardPatternTest, GetIndexReturnsCorrectIndexWithMockSiteAndParamType )
{
	// for this pattern the lattice size is irrelevant (since lattice index is running slowest and no parity splitting)
	SiteTypeMock<> mySite( 5 );
	const int myParamSize = 18;
	const int myParamIndex = 4;
	const int nDim = 4;
	const int myMu = 1;
	int expect = (mySite.index*nDim+myMu)*myParamSize+myParamIndex;

	int result = StandardPattern<SiteTypeMock<>, ParamTypeMock<myParamSize> >::getIndex( mySite, myMu, myParamIndex );

	ASSERT_EQ( expect, result );
}

TEST( StandardPatternTest, GetIndexIsCompatibleToGetStandardIndex )
{
	const int siteIndex = 5;
	const int nDim = 4;
	SiteTypeMock<nDim> mySite( siteIndex );
	const int myParamSize = 18;
	const int myParamIndex = 4;
	const int myMu = 1;
	int expect = StandardPattern<SiteTypeMock<nDim>, ParamTypeMock<myParamSize> >::getIndex( mySite, myMu, myParamIndex );

	int result = StandardPattern<SiteTypeMock<nDim>, ParamTypeMock<myParamSize> >::getStandardIndex( siteIndex, myMu, myParamIndex );

	ASSERT_EQ( expect, result );
}

TEST( StandardPatternTest, PatternHasDefinedTypeSITETYPE )
{
	typename StandardPattern<SiteTypeMock<4>, ParamTypeMock<2> >::SITETYPE s(5);

	ASSERT_EQ( 5, s.getIndex() );
}
