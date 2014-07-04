


#include "gmock/gmock.h"
#include "LinkToGaugefieldConverter.h"
#include "LocalLink.h"
#include "../common/culgt_typedefs.h"
#include "parameterization_types/SUNComplexFull.h"
#include "../cuLGT1legacy/Complex.hxx"
using namespace culgt;
using namespace ::testing;


void ASSERT_NEAR_FLOAT4( float4 expected, float4 actual, double precision = 1E-7)
{
	ASSERT_NEAR( expected.x, actual.x, precision );
	ASSERT_NEAR( expected.y, actual.y, precision );
	ASSERT_NEAR( expected.z, actual.z, precision );
	ASSERT_NEAR( expected.w, actual.w, precision );
}
void ASSERT_NEAR_DOUBLE4( double4 expected, double4 actual, double precision = 1E-14 )
{
	ASSERT_NEAR( expected.x, actual.x, precision );
	ASSERT_NEAR( expected.y, actual.y, precision );
	ASSERT_NEAR( expected.z, actual.z, precision );
	ASSERT_NEAR( expected.w, actual.w, precision );
}

void ASSERT_NEAR_EQ_ARRAY( double* expected, double* actual, int size, double precision = 1E-14 )
{
	for( int i = 0; i < size; i++ )
	{
		ASSERT_NEAR( expected[i], actual[i], precision );
	}
}

class LinkToGaugefieldConverterSU2: public Test
{
public:
	LocalLink<SU2Vector4<double> > link;
	Real4<double>::VECTORTYPE linkData;
	const double someValue;

	LinkToGaugefieldConverterSU2() : someValue( .5 )
	{
		linkData.x = someValue;
		linkData.y = someValue;
		linkData.z = someValue;
		linkData.w = someValue;
		link.set(0, linkData );
	}
};

TEST_F( LinkToGaugefieldConverterSU2, ExtractsLinearGaugefield )
{
	double A[3];

	LinkToGaugefieldConverter<2,gaugefieldtype::LINEAR>::convert( A, link );

	ASSERT_DOUBLE_EQ( 2.*someValue, A[0] );
}

TEST_F( LinkToGaugefieldConverterSU2, ConvertsLinearGaugefieldBackToLink )
{
	double A[3];
	LinkToGaugefieldConverter<2,gaugefieldtype::LINEAR>::convert( A, link );

	LocalLink<SU2Vector4<double> > link2;

	LinkToGaugefieldConverter<2,gaugefieldtype::LINEAR>::convert( link2, A );

	ASSERT_TRUE( link2 == link ); // TODO this should be a ASSERT_FLOAT style comparison.
}

TEST_F( LinkToGaugefieldConverterSU2, ConvertsLogarithmicGaugefieldBackToLink )
{
	double A[3];
	LinkToGaugefieldConverter<2,gaugefieldtype::LOGARITHMIC>::convert( A, link );

	LocalLink<SU2Vector4<double> > link2;
	link2.zero();

	LinkToGaugefieldConverter<2,gaugefieldtype::LOGARITHMIC>::convert( link2, A );

	ASSERT_NEAR_DOUBLE4( link.get(0), link2.get(0) );
}


class LinkToGaugefieldConverterSU3: public Test
{
public:
	LocalLink<SUNComplexFull<3,double> > link;
	const double someValue;
	double Alin[8];
	double Alog[8];
	double A[8];

	typedef Complex<double> MyComplex;

	LinkToGaugefieldConverterSU3() : someValue( .3124 )
	{
		// IMPORTANT: the exact values are calculated with Matlab but truncated, that's why we don't reach full double precisions

		link << MyComplex(0.027288592117974 , 0.340040894989286), 	MyComplex(-0.076215048176904 , 0.636601268202080),  MyComplex(0.516066183900287 , 0.454129175526738),
				MyComplex(-0.849220194098520 , 0.308910080079893),  MyComplex(0.127532894766798, - 0.211665307965350),  MyComplex(0.291452572369848, - 0.193360235249234),
				MyComplex(-0.061058506982911, - 0.251592936229194), MyComplex(-0.652242622656688 , 0.320066885355475),  MyComplex(0.207384244638884, - 0.601748600604039);
		Alin[0] = 0.945511348281973;
		Alin[1] = 0.773005145921617;
		Alin[2] = 0.551706202954636;
		Alin[3] = 0.202536239297544;
		Alin[4] = 0.577124690883198;
		Alin[5] = 0.126706650106241;
		Alin[6] = 0.943695195026536;
		Alin[7] = 0.768957112812092;

		Alog[0] = 1.837681126110123;
		Alog[1] = 1.523400238601420;
		Alog[2] = 1.093818132657546;
		Alog[3] = 0.425629843562787;
		Alog[4] = 1.114794301271606;
		Alog[5] = 0.234672350073825;
		Alog[6] = 1.904512632090502;
		Alog[7] = 1.524034585031251;
	}
};

TEST_F( LinkToGaugefieldConverterSU3, ExtractsLinearGaugefield )
{
	LinkToGaugefieldConverter<3,gaugefieldtype::LINEAR>::convert( A, link );

	ASSERT_NEAR_EQ_ARRAY( Alin, A, 8 );
}

TEST_F( LinkToGaugefieldConverterSU3, ExtractsLogarithmicGaugefield )
{
	LinkToGaugefieldConverter<3,gaugefieldtype::LOGARITHMIC>::convert( A, link );

	ASSERT_NEAR_EQ_ARRAY( Alog, A, 8 );
}





