#include "gmock/gmock.h"
#include "cudatest/cuda_gtest_plugin.h"
#include "cuLGT1legacy/Complex.hxx"
using namespace ::testing;
using namespace culgt;


CUDA_TEST(AComplex, SizeOfComplexDoubleIs16 )
{
	int expect = 16;

	ASSERT_EQ( expect, (int)sizeof( Complex<double> ) );
}

CUDA_TEST(AComplex, SizeOfComplexFloatIs8 )
{
	int expect = 8;

	ASSERT_EQ( expect, (int)sizeof( Complex<float> ) );
}

CUDA_TEST(AComplex, MultByI )
{
	Complex<float> temp(1.,2.);

	temp *= Complex<float>::I();

	ASSERT_FLOAT_EQ( (float)-2., temp.x );
	ASSERT_FLOAT_EQ( (float)1., temp.y );
}

CUDA_TEST(AComplex, Sqrt )
{
	Complex<float> temp(1.,2.);

	Complex<float>temp2;
	temp2 = sqrt(temp);

	temp2 *= temp2;

	ASSERT_FLOAT_EQ( temp.x, temp2.x );
	ASSERT_FLOAT_EQ( temp.y, temp2.y );
}

CUDA_TEST(AComplex, Log )
{
	Complex<float> temp(1.,2.);
	Complex<float> result(  8.0471897e-01, 1.1071488e+00 );

	temp = log(temp);

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

CUDA_TEST(AComplex, Acos )
{
	Complex<float> temp(1.,2.);

	Complex<float> result( 1.1437178e+00, - 1.5285709e+00);

	temp = acos(temp);

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

CUDA_TEST(AComplex, Asin )
{
	Complex<float> temp(1.,2.);

	Complex<float> result( 4.2707860e-01, 1.5285709e+00);

	temp = asin(temp);

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

CUDA_TEST(AComplex, CosOfAcos )
{
	Complex<float> result(1.,2.);

	Complex<float> temp( result.x, result.y);
	temp = cos(acos(temp));

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

CUDA_TEST(AComplex, AsinForProblematicCaseWhereLogDefOfAcosFails )
{
	Complex<float> temp(-27309.294921875000000,-22527.775390625000000);
	Complex<float> result(-0.8810484, -11.1676693);

	temp = asin( temp );

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

/**
 * fails in host code (numerical inaccuracy)
 */
CUDA_TEST(AComplex, DISABLED_AsinForProblematicCaseWhereLogDefOfAcosFails2 )
{
	Complex<float> temp(-1131.254760742187500,4343.823730468750000);
	Complex<float> result( -2.5476924e-01, 9.1024685e+00);


	temp = asin( temp );

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}

CUDA_TEST(AComplex, AsinForProblematicCaseWhereLogDefOfAcosFails3 )
{
	Complex<float> temp(-5222.564453125000000,-9092.003906250000000);
	Complex<float> result( 2.092189281608619e+00,9.950868543459121e+00);


	temp = acos( temp );

	ASSERT_FLOAT_EQ( result.x, temp.x );
	ASSERT_FLOAT_EQ( result.y, temp.y );
}
