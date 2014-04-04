

#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"
#include "GPUPatternTimesliceParityPriority.h"
#include "../LatticeDimension.h"

using namespace culgt;
using namespace ::testing;



TEST( GPUPatternTimesliceParityPriorityTest, GetIndexReturnsCorrectIndex )
{
	const int Ndim = 4;
	LatticeDimension<Ndim> dim( 4, 4, 4, 4 );
	int latSizeTimeslice = dim.getSizeTimeslice();
	const int siteIndexTimeslice = 35%(latSizeTimeslice/2);
	const int siteIndexTimesliceInParity = siteIndexTimeslice%(latSizeTimeslice/2);
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const int t = 2;
	int siteIndex = siteIndexTimeslice+t*latSizeTimeslice;

	int parity = (siteIndexTimeslice/(latSizeTimeslice/2)); // = 1

	SiteTypeMock<Ndim,TIMESLICE_SPLIT> mySite(dim.getSize(),latSizeTimeslice,siteIndex,siteIndexTimeslice);

	int expect = (((t*2+parity)*Ndim+myMu)*paramTypeSize+paramIndex)*(latSizeTimeslice/2)+siteIndexTimesliceInParity;

	int result = GPUPatternTimesliceParityPriority<SiteTypeMock<Ndim,TIMESLICE_SPLIT>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	ASSERT_EQ( expect, result );
}
