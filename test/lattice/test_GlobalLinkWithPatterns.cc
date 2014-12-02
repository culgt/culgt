#include "gmock/gmock.h"

#include "lattice/GlobalLink.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/configuration_patterns/StandardPattern.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/SU3Real12.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"

using namespace std;

using namespace culgt;
using namespace ::testing;

typedef SiteIndex<4,NO_SPLIT> SiteIndex4;
typedef GlobalLink<StandardPattern<SiteIndex4,SUNRealFull<3,float> > > GlobalLinkStandardSU3Real18Index;
typedef GlobalLink<StandardPattern<SiteIndex4,SU3Real12<float> > > GlobalLinkStandardSU3Real12Index;


void ASSERT_EQ_ZERO( GlobalLinkStandardSU3Real18Index link )
{
	for( int i = 0; i < GlobalLinkStandardSU3Real18Index::CONFIGURATIONPATTERN::PARAMTYPE::SIZE; i++ )
	{
		ASSERT_FLOAT_EQ( 0., link.get(i) );
	}
}

class AGlobalLinkWithPattern: public Test
{
public:
	static const int size[4];
	SiteIndex4 s;
	SiteIndex4 s2;
	static const int mu = 0;
	static const int parameterIndex = 0;
	float U[2*2*2*2*4*18];
	float someValue;

	AGlobalLinkWithPattern() : s(size), s2(size), someValue(1.234)
	{
	}
};

const int AGlobalLinkWithPattern::size[4] = {2,2,2,2};


TEST_F(AGlobalLinkWithPattern, TypedefSITETYPEworks )
{
	GlobalLinkStandardSU3Real18Index::CONFIGURATIONPATTERN::SITETYPE s(size);
	s.setIndex( 1 );
	ASSERT_EQ( 1, s.getIndex() );
}

TEST_F(AGlobalLinkWithPattern, SetGetValueWithSameLinkAndParameterIndexWorks )
{
	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );

	link.set( parameterIndex, someValue );

	ASSERT_FLOAT_EQ( someValue, link.get( parameterIndex) );
}

TEST_F(AGlobalLinkWithPattern, ZeroSetsAllElementsOfLinkToZero )
{
	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );

	link.set( 1, someValue );
	link.zero();

	ASSERT_EQ_ZERO( link );
}

TEST_F(AGlobalLinkWithPattern, OperatorAssignWithSameTypeDoesNotChangeSite )
{
	int link2Index = 3;

	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );
	s2.setIndex( link2Index );
	GlobalLinkStandardSU3Real18Index link2( U, s2, mu );

	link2 = link;

	ASSERT_EQ( link2Index, link2.site.getIndex() );
}

TEST_F(AGlobalLinkWithPattern, OperatorAssignWithSameTypeCopiesData )
{
	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );
	s2.setIndex( 3 );
	GlobalLinkStandardSU3Real18Index link2( U, s2, mu );

	link.set( 1, someValue );

	link2.zero();
	link2 = link;

	ASSERT_FLOAT_EQ( someValue, link2.get( 1 ) );
}

class AGlobalLinkWithSU3Real12AndSU3Real18: public Test
{
public:
	static const int size[4];
	SiteIndex4 s;
	SiteIndex4 s2;
	static const int mu = 0;
	static const int parameterIndex = 0;
	float U[2*2*2*2*4*18];
	float someValue;
	GlobalLinkStandardSU3Real18Index* link18;
	GlobalLinkStandardSU3Real12Index* link12;

	AGlobalLinkWithSU3Real12AndSU3Real18() : s(size), s2(size), someValue(1.234)
	{
		s.setLatticeIndex( 0 );
		link18 = new GlobalLinkStandardSU3Real18Index( U, s, mu );
		link18->zero();
		s2.setLatticeIndex( 1 );
		link12 = new GlobalLinkStandardSU3Real12Index( U, s2, mu );
		link12->zero();
	}
};
const int AGlobalLinkWithSU3Real12AndSU3Real18::size[4] = {2,2,2,2};

TEST_F( AGlobalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesFromReal18ToReal12 )
{
	link18->set(1, 2.42);
	*link12 = *link18;
	ASSERT_FLOAT_EQ( 2.42, link12->get(1) );
}

TEST_F( AGlobalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLine )
{
	// define a matrix that has (1 0 0) in first and (0 1 0) in second line: expect (0 0 1) in third line
	link12->set(0, 1.);
	link12->set(8, 1.);

	*link18 = *link12;

	ASSERT_FLOAT_EQ( 1., link18->get(16) );
}
