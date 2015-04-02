#include "gmock/gmock.h"
#include <iostream>
#include <iomanip>
#include "math/Complex.h"
#include "eigensolver/matrix_vector.h"
#include "eigensolver/Cardano3x3.h"

using namespace culgt;
using namespace testing;
using std::setprecision;

//TEST( Cardano3x3Double, ComputesEigenvaluesForDiagonalMatrix )
//{
//	const double diagonal[3] = {.5,1.,2.};
//	Matrix<double,3> matrix( diagonal[0],diagonal[1], diagonal[2] );
//	Cardano3x3<double> eigen( matrix );
//
//	ASSERT_DOUBLE_EQ( diagonal[0], eigen.getEigenvalue( 0 ) );
//	ASSERT_DOUBLE_EQ( diagonal[1], eigen.getEigenvalue( 1 ) );
//	ASSERT_DOUBLE_EQ( diagonal[2], eigen.getEigenvalue( 2 ) );
//}
//
//TEST( Cardano3x3Double, ComputesEigenvaluesForEasyMatrix )
//{
//	const double diagonal[3] = {.5,.1,.2};
//	Matrix<double,3> matrix( diagonal[0],diagonal[1], diagonal[2] );
//	matrix(2,0) = .5;
//	Cardano3x3<double> eigen( matrix );
//
//	ASSERT_DOUBLE_EQ( diagonal[0], eigen.getEigenvalue( 0 ) );
//	ASSERT_DOUBLE_EQ( diagonal[1], eigen.getEigenvalue( 1 ) );
//	ASSERT_DOUBLE_EQ( diagonal[2], eigen.getEigenvalue( 2 ) );
//}
//
//TEST( Cardano3x3Double, ComputesEigenvaluesForFullMatrix )
//{
//	double values[9] = {0.5648000,   0.9337000 ,  0.5146000,
//			   0.0044000 ,  0.7516000 ,  0.0349000,
//			   0.9919000,   0.7123000 ,  0.7573000 };
//	double eigenvals[3] = {    1.438722545459362e+00,
//		    -4.629879539880649e-02,
//		     6.812762499394441e-01};
//
//	Matrix<double,3> matrix( values );
//
//	Cardano3x3<double> eigen( matrix );
//
//	ASSERT_FLOAT_EQ( eigenvals[0], eigen.getEigenvalue( 0 ) );
//	ASSERT_FLOAT_EQ( eigenvals[1], eigen.getEigenvalue( 1 ) );
//	ASSERT_FLOAT_EQ( eigenvals[2], eigen.getEigenvalue( 2 ) );
//}


class Cardano3x3Complex: public Test
{
public:
	typedef Complex<double> MyComplex;
	MyComplex values[9];
	MyComplex eigenvals[3];

	Matrix<MyComplex,3> matrix;

	Cardano3x3Complex()
	{
		matrix(0,0) = MyComplex(0.0373 , 0.7375);
		matrix(0,1) = MyComplex(0.1872 , 0.9202);
		matrix(0,2) = MyComplex(0.8597 , 0.1376);
		matrix(1,0) = MyComplex(0.8900 , 0.6943);
		matrix(1,1) = MyComplex(0.6674 , 0.3166);
		matrix(1,2) = MyComplex(0.7268 , 0.0941);
		matrix(2,0) = MyComplex(0.0872 , 0.9044);
		matrix(2,1) = MyComplex(0.3103 , 0.4665);
		matrix(2,2) = MyComplex(0.0552 , 0.5392);

		eigenvals[0] = MyComplex( 1.311611535467893e+00, + 1.749685786245779e+00 );
		eigenvals[1] = MyComplex( -4.479438961528297e-01, - 1.863852435495146e-01 );
		eigenvals[2] = MyComplex( -1.037676393150646e-01, + 2.999945730373376e-02 );

	}
};

/*
 * We do the calculation in double, but compare with FLOAT_EQ because of rounding errors.
 * If you don't know what to do else: fix this...
 */
TEST_F( Cardano3x3Complex, ComputesEigenvaluesForFullMatrix )
{
	Cardano3x3<MyComplex> eigen( matrix );

	ASSERT_FLOAT_EQ( eigenvals[0].x, eigen.getEigenvalue( 0 ).x );
	ASSERT_FLOAT_EQ( eigenvals[0].y, eigen.getEigenvalue( 0 ).y );
	ASSERT_FLOAT_EQ( eigenvals[1].x, eigen.getEigenvalue( 1 ).x );
	ASSERT_FLOAT_EQ( eigenvals[1].y, eigen.getEigenvalue( 1 ).y );
	ASSERT_FLOAT_EQ( eigenvals[2].x, eigen.getEigenvalue( 2 ).x );
	ASSERT_FLOAT_EQ( eigenvals[2].y, eigen.getEigenvalue( 2 ).y );
}

TEST_F( Cardano3x3Complex, ComputesEigenvectorsForFullMatrix )
{
	Cardano3x3<MyComplex> eigen( matrix );


	Vector<MyComplex,3> result;

	for( int i = 0; i < 3; i++ )
	{
		result = eigen.getEigenvector(i);

		Vector<MyComplex,3> check;
		check = (matrix*result);
		check -= result*eigen.getEigenvalue(i);

		ASSERT_NEAR( 0., check( 0 ).x, 1e-14 );
		ASSERT_NEAR( 0., check( 0 ).y, 1e-14 );
		ASSERT_NEAR( 0., check( 1 ).x, 1e-14 );
		ASSERT_NEAR( 0., check( 1 ).y, 1e-14 );
		ASSERT_NEAR( 0., check( 2 ).x, 1e-14 );
		ASSERT_NEAR( 0., check( 2 ).y, 1e-14 );
	}
}



class Cardano3x3Unitary: public Test
{
public:
	typedef Complex<double> MyComplex;
	MyComplex values[9];
	MyComplex eigenvals[3];

	Matrix<MyComplex,3> matrix;

	Cardano3x3Unitary()
	{
		matrix(0,0) = MyComplex(2.728859211797367e-02, + 3.400408949892865e-01);
		matrix(0,1) = MyComplex(-7.621504817690354e-02, + 6.366012682020801e-01);
		matrix(0,2) = MyComplex(5.160661839002871e-01, + 4.541291755267379e-01);
		matrix(1,0) = MyComplex(-8.492201940985201e-01, + 3.089100800798935e-01);
		matrix(1,1) = MyComplex(1.275328947667976e-01, - 2.116653079653496e-01);
		matrix(1,2) = MyComplex(2.914525723698476e-01, - 1.933602352492340e-01);
		matrix(2,0) = MyComplex(-6.105850698291067e-02, - 2.515929362291943e-01);
		matrix(2,1) = MyComplex(-6.522426226566884e-01, + 3.200668853554753e-01);
		matrix(2,2) = MyComplex(2.073842446388842e-01, - 6.017486006040392e-01);
	}
};

TEST_F( Cardano3x3Unitary, ComputesLogOfMatrix )
{
	Cardano3x3<MyComplex> eigen( matrix );

	MyComplex expect00(  6.036837696399289e-16, 9.868599552898188e-01 );
	MyComplex expect21( -9.522563160452506e-01, 1.173361750369123e-01 );

	eigen.logAssumeUnitary();

	for( int i = 0; i < 3; i++ )
	{
		ASSERT_NEAR( expect00.x, matrix(0,0).x, 1e-14 );
		ASSERT_NEAR( expect00.y, matrix(0,0).y, 1e-14);
		ASSERT_NEAR( expect21.x, matrix(2,1).x, 1e-14 );
		ASSERT_NEAR( expect21.y, matrix(2,1).y, 1e-14 );
	}
}
