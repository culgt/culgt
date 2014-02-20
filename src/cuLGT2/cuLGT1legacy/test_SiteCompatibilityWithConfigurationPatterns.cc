/**
 * test_SiteIndex.cc
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */



#include "gmock/gmock.h"
#include "../lattice/configuration_patterns/StandardPattern.h"
#include "../lattice/configuration_patterns/testhelper_PatternMocks.h"
#include "SiteIndex.hxx"
#include "SiteCoord.hxx"
#include "../lattice/configuration_patterns/testhelper_PatternMocks.h"

using namespace ::testing;
using namespace culgt;


class StandardPatternCompatible: public Test
{
public:
	static const int nDim = 4;

	const int latIndex = 5;
	int size[nDim];
	SiteIndex<nDim,NO_SPLIT> mySite;

	StandardPatternCompatible(): size{4,4,4,4}, mySite(size){};
};

TEST_F(StandardPatternCompatible, WithSiteIndexNoSplit )
{
	const int myParamIndex = 4;
	mySite.setLatticeIndex( latIndex );

	int result = StandardPattern<SiteIndex<nDim,NO_SPLIT>, ParamTypeMock<1> >::getIndex( mySite, 0, myParamIndex );

	ASSERT_EQ( latIndex*4+myParamIndex, result );
}
