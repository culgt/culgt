

#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"
#include "GPUPatternTimeslice.h"
#include "../LatticeDimension.h"

using namespace culgt;
using namespace ::testing;



TEST( GPUPatternTimesliceTest, GetIndexReturnsCorrectIndex )
{
	const int Ndim = 4;
	LatticeDimension<Ndim> dim( 4, 4, 4, 4 );
	int latSizeTimeslice = dim.getSizeTimeslice();
	const int siteIndexTimeslice = 10;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const int t = 2;
	int siteIndex = siteIndexTimeslice+t*latSizeTimeslice;

	SiteTypeMock<Ndim> mySite(dim.getSize(),latSizeTimeslice,siteIndex,siteIndexTimeslice);

	int expect = ((t*Ndim+myMu)*paramTypeSize+paramIndex)*latSizeTimeslice+siteIndexTimeslice;

	int result = GPUPatternTimeslice<SiteTypeMock<Ndim>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPatternTimesliceTest, ConvertToStandardIndexIsCompatibleToStandardPattern )
{
	const int Ndim = 4;
	LatticeDimension<Ndim> dim( 4, 4, 4, 4 );
	int latSizeTimeslice = dim.getSizeTimeslice();
	const int siteIndexTimeslice = 10;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const int t = 2;
	int siteIndex = siteIndexTimeslice+t*latSizeTimeslice;

	SiteTypeMock<Ndim> mySite(dim.getSize(),latSizeTimeslice,siteIndex,siteIndexTimeslice);

	int expect = StandardPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	int gpuPatternIndex = GPUPatternTimeslice<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );
	int result = GPUPatternTimeslice<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::convertToStandardIndex( gpuPatternIndex, dim );

	ASSERT_EQ( expect, result );
}
