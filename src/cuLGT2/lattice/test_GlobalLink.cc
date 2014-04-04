
#include "GlobalLink.h"
#include "gmock/gmock.h"

using namespace culgt;
using namespace ::testing;


class ParameterTypeStub
{
public:
	typedef float TYPE;
};

class SiteStub
{
public:
	static const int Ndim = 4;
	int index;
	SiteStub() : index(0) {}
	SiteStub( int i ) : index(i) {}
};

template<int N=0> class AccessPatternStub
{
public:
	typedef SiteStub SITETYPE;
	typedef ParameterTypeStub PARAMTYPE;

	static int getIndex( SiteStub site, int mu, int patternIndex )
	{
		return N;
	}
};

class AGlobalLink: public Test
{
public:
	SiteStub s;
	static const int mu = 0;
	static const int parameterIndex = 0;
	float U[100];
};

TEST_F(AGlobalLink, SetGetValueWithSameLinkAndParameterIndexWorks )
{
	GlobalLink<AccessPatternStub<> > link( U, s, mu );
	link.set( parameterIndex, 1.234 );

	ASSERT_FLOAT_EQ( 1.234, link.get( parameterIndex ) );
}

TEST_F(AGlobalLink, GetsValueFromTheCorrectPosition )
{
	U[5] = 1.234;
	GlobalLink<AccessPatternStub<5> > link( U, s, mu );

	ASSERT_FLOAT_EQ( 1.234, link.get(parameterIndex) );
}

TEST_F(AGlobalLink, SetsValueToTheCorrectPosition )
{
	U[5] = -1.;
	GlobalLink<AccessPatternStub<5> > link( U, s, mu );
	link.set( parameterIndex, 1.234 );

	ASSERT_FLOAT_EQ( 1.234, U[5] );
}

TEST_F(AGlobalLink, SiteIsKeptAsCopy )
{
	SiteStub site;
	site.index = 5;
	GlobalLink<AccessPatternStub<5> > link( U, site, mu );
	site.index = 3;

	ASSERT_EQ( 5, link.site.index );
}

TEST(GlobalLinkTemplateParameterTest, SetsSITETYPE )
{
	typename GlobalLink<AccessPatternStub<5> >::CONFIGURATIONPATTERN::SITETYPE s(1);
	ASSERT_EQ( 1, s.index );
}

