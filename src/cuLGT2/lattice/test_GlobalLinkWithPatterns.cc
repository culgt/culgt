#include "gmock/gmock.h"

#include "GlobalLink.h"
#include "../cuLGT1legacy/SiteIndex.hxx"
#include "configuration_patterns/StandardPattern.h"
#include "su3/SU3Real18.h"

using namespace std;

using namespace culgt;
using namespace ::testing;


class AGlobalLinkWithPattern: public Test
{
public:
	typedef SiteIndex<4,NO_SPLIT> SiteIndex4;
	static const int size[4];
	SiteIndex4 s;
	static const int mu = 0;
	static const int parameterIndex = 0;


	AGlobalLinkWithPattern() : s(size)
	{
	}
};

const int AGlobalLinkWithPattern::size[4] = {2,2,2,2};

TEST_F(AGlobalLinkWithPattern, SetsSITETYPE )
{
	GlobalLink<StandardPattern<SiteIndex4,SU3Real18<float> > >::CONFIGURATIONPATTERN::SITETYPE s(size);
	s.setIndex( 1 );
	ASSERT_EQ( 1, s.getIndex() );
}

TEST_F(AGlobalLinkWithPattern, SetGetValueWithSameLinkAndParameterIndexWorks )
{
	float U[2*2*2*2*4*18];
	s.setLatticeIndex( 1 );
	GlobalLink<StandardPattern<SiteIndex4,SU3Real18<float> > > link( U, s, mu );

	link.set( parameterIndex, 1.234 );

	ASSERT_FLOAT_EQ( 1.234, link.get( parameterIndex) );
}

