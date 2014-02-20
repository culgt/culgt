#include "gmock/gmock.h"
#include "LocalLink.h"
#include "su3/SU3Real12.h"
#include "su3/SU3Real18.h"
#include "su3/ParameterizationMediatorSU3_Real12_Real18.h"

using namespace culgt;
using namespace ::testing;


class ALocalLinkWithSU3Real12AndSU3Real18: public Test
{
public:
	LocalLink<SU3Real18<float> > link18;
	LocalLink<SU3Real12<float> > link12_1;
	LocalLink<SU3Real12<float> > link12_2;

	void SetUp()
	{
		link18.zero();
		link12_1.zero();
		link12_2.zero();
	}
};

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, GetReturnsPreviouslySetValue )
{
	link18.set(1, 2.42);
	ASSERT_FLOAT_EQ( 2.42, link18.get(1) );
}

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesForSameParamType )
{
	link12_1.set(1, 2.42);
	link12_2 = link12_1;
	ASSERT_FLOAT_EQ( 2.42, link12_2.get(1) );
}

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, OperatorAssignDoesNotUseReference )
{
	link12_1.set(1, 2.42);
	link12_2 = link12_1;
	link12_2.set(1, 1.3 );
	EXPECT_FLOAT_EQ( 1.3, link12_2.get(1) );
	ASSERT_FLOAT_EQ( 2.42, link12_1.get(1) );
}

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesFromReal18ToReal12 )
{
	link18.set(1, 2.42);
	link12_1 = link18;
	ASSERT_FLOAT_EQ( 2.42, link12_1.get(1) );
}

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesFromReal12ToReal18 )
{
	link12_1.set(1, 2.42);
	link18 = link12_1;
	ASSERT_FLOAT_EQ( 2.42, link18.get(1) );
}

TEST_F( ALocalLinkWithSU3Real12AndSU3Real18, OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLine )
{
	// define a matrix that has (1 0 0) in first and (0 1 0) in second line: expect (0 0 1) in third line
	link12_1.set(0, 1.);
	link12_1.set(8, 1.);
	link18 = link12_1;
	ASSERT_FLOAT_EQ( 1., link18.get(16) );
}
