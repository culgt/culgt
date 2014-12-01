

#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"
#include "GPUPattern.h"

using namespace culgt;
using namespace ::testing;

TEST( GPUPatternTest, GetIndexReturnsSiteIndexForParamIndexZeroAndMuZero )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	SiteTypeMock<> mySite(latSize, siteIndex);

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<0> >::getIndex( mySite, 0, 0 );

	ASSERT_EQ( mySite.index, result );
}

TEST( GPUPatternTest, GetIndexReturnsCorrectIndexWithMuZero )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramIndex = 4;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = paramIndex*latSize+siteIndex;

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, 0, paramIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPatternTest, GetIndexReturnsCorrectIndex )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = (myMu*paramTypeSize+paramIndex)*latSize+siteIndex;

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPatternTest, ConvertToStandardIndexIsCompatibleToStandardPattern )
{
	LatticeDimension<4> dim( 1, 10, 10, 10 );
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = StandardPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	int gpuPatternIndex = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, myMu, paramIndex );
	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::convertToStandardIndex( gpuPatternIndex, dim );

	ASSERT_EQ( expect, result );
}
