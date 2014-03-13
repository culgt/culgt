#include "gmock/gmock.h"
#include "LocalLink.h"
#include "parameterization_types/SU3Real12.h"
#include "parameterization_types/SUNRealFull.h"
#include "parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"

using namespace culgt;
using namespace ::testing;


class ALocalLinkWithSU3Real12AndSU3Real18: public Test
{
public:
	LocalLink<SUNRealFull<3,float> > link18;
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

class ASpecialLocalLinkWithSU3Real18: public Test
{
public:
	LocalLink<SUNRealFull<3,float> > link18;

	const float redet = 0.585224903024000;
	const float imdet = 0.185504451600000;
	const float retrace = 2.3660;
	const float imtrace = 2.3844;

	void SetUp()
	{
		/*
		 * set a 3x3 matrix A with
		 * det(A) =  0.585224903024000 - 0.185504451600000i
		 * trace(A) = 2.3660 + 2.3844i
		 */

		link18.set(0, 0.9649);
		link18.set(1, 0.7922);
		link18.set(2, 0.9572);
		link18.set(3, 0.0357);
		link18.set(4, 0.1419);
		link18.set(5, 0.6787);
		link18.set(6, 0.1576);
		link18.set(7, 0.9595);
		link18.set(8, 0.4854);
		link18.set(9, 0.8491);
		link18.set(10, 0.4218);
		link18.set(11, 0.7577);
		link18.set(12, 0.9706);
		link18.set(13, 0.6557);
		link18.set(14, 0.8003);
		link18.set(15, 0.9340);
		link18.set(16, 0.9157);
		link18.set(17, 0.7431);
	}
};

TEST_F( ASpecialLocalLinkWithSU3Real18, IdentityWorks )
{
	link18.identity();

	ASSERT_FLOAT_EQ( 1., link18.get(0) );
}

TEST_F( ASpecialLocalLinkWithSU3Real18, ReDetWorks )
{
	float result = link18.reDet();

	ASSERT_FLOAT_EQ( redet, result );
}

TEST_F( ASpecialLocalLinkWithSU3Real18, ReTraceWorks )
{
	float result = link18.reTrace();

	ASSERT_FLOAT_EQ( retrace, result );
}





