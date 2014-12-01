#include "gmock/gmock.h"
#include "../configuration_patterns/testhelper_PatternMocks.h"
#include "DiracAdjointColorFieldStructure.h"

using namespace culgt;
using namespace ::testing;

TEST( ADiracAdjointColorFieldStructure, ConstructorWithSiteIndexOnlyEqualsFieldIndex )
{
	SiteTypeMock<> site(32);
	DiracAdjointColorFieldStructure<SiteTypeMock<>, 3> fieldStructure( site );

	ASSERT_EQ( site.getIndex() , fieldStructure.getIndex() );
}

TEST( ADiracAdjointColorFieldStructure, ChangingSiteIndexWorks )
{
	SiteTypeMock<> site(0);
	const int latIndex = 344;
	DiracAdjointColorFieldStructure<SiteTypeMock<>, 3> fieldStructure( site );

	fieldStructure.getSite().setIndex( latIndex );

	ASSERT_EQ( latIndex , fieldStructure.getIndex() );
}

TEST( ADiracAdjointColorFieldStructure, FullConstructorWorks )
{
	int latSize = 10;
	const int NDIM = 4;
	const int NC = 3;
	const int latIndex = 5;
	const int a = 1;
	const int mu = 2;

	SiteTypeMock<NDIM> site(latSize,latIndex);
	DiracAdjointColorFieldStructure<SiteTypeMock<NDIM>, NC> fieldStructure( mu, a, site );

	ASSERT_EQ( (mu*(NC*NC-1)+a)*latSize+latIndex , fieldStructure.getIndex() );
}
