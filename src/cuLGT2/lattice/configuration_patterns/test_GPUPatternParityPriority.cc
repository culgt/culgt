
#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"
#include "GPUPatternParityPriority.h"

using namespace culgt;
using namespace ::testing;

TEST( GPUPatternParityPriorityTest, GetIndexReturnsSiteIndexForIndexInLowerHalfOfLattice )
{
	const int paramTypeSize = 10;
	const int latSize = 200;
	const int siteIndex = 10;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int result = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, 0, 0 );

	ASSERT_EQ( mySite.index, result );
}

TEST( GPUPatternParityPriorityTest, GetIndexReturnsSiteIndexForIndexInUpperHalfOfLattice )
{
	const int paramTypeSize = 10;
	const int paramTypeIndex = 1;
	const int latSize = 200;
	const int siteIndex = 110;
	SiteTypeMock<1> mySite(latSize,siteIndex);

	//expect:
	int parity = (siteIndex/(latSize/2)); // = 1
	int siteIndexInParity = (siteIndex%(latSize/2));

	int expect = (parity*paramTypeSize+paramTypeIndex)*(latSize/2)+siteIndexInParity;

	int result = GPUPatternParityPriority<SiteTypeMock<1>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, 0, paramTypeIndex );

	ASSERT_EQ( expect, result );
}


TEST( GPUPatternParityPriorityTest, GetIndexReturnsSiteIndexForIndexInUpperHalfOfLatticeWithMu )
{
	const int paramTypeSize = 10;
	const int paramTypeIndex = 1;
	const int latSize = 200;
	const int siteIndex = 110;
	const int nDim = 4;
	const int myMu = 2;
	SiteTypeMock<nDim> mySite(latSize,siteIndex);

	//expect:
	int parity = (siteIndex/(latSize/2)); // = 1
	int siteIndexInParity = (siteIndex%(latSize/2));

	int expect = ((parity*nDim+myMu)*paramTypeSize+paramTypeIndex)*(latSize/2)+siteIndexInParity;

	int result = GPUPatternParityPriority<SiteTypeMock<nDim>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramTypeIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPatternParityPriorityTest, ConvertToStandardIndexIsCompatibleToStandardPattern )
{
	const int latSize = 200;
	const int siteIndex = 110;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = StandardPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	int gpuPatternIndex = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, myMu, paramIndex );
	int result = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<18> >::convertToStandardIndex( gpuPatternIndex, mySite.getSize() );

	ASSERT_EQ( expect, result );
}
