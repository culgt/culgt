#include "eigensolver/Matrix.h"
#include "gmock/gmock.h"

using namespace culgt;

TEST( Matrix3x3Double, ElementAccessOperatorWorks )
{
	const double someValue = 1.2423;
	Matrix<double,3> matrix;

	matrix( 2, 1 ) = someValue;

	ASSERT_DOUBLE_EQ( someValue, matrix(2,1) );
}

TEST( Matrix3x3Double, Trace )
{
	const double someValue = 1.2423;
	Matrix<double,3> matrix( someValue, someValue, someValue );

	ASSERT_DOUBLE_EQ( 3.*someValue, matrix.trace() );
}

TEST( Matrix3x3Double, Subtract )
{
	const double someValue = 1.2423;
	const double someValue2 = 4.2332;
	Matrix<double,3> matrix( someValue, someValue, someValue );

	Matrix<double,3> matrix2;

	matrix2 = matrix-someValue2;

	ASSERT_DOUBLE_EQ( someValue-someValue2, matrix2(0,0) );
}

TEST( Matrix3x3Double, Divide )
{
	const double someValue = 1.2423;
	const double someValue2 = 4.2332;
	Matrix<double,3> matrix( someValue, someValue, someValue );

	Matrix<double,3> matrix2;

	matrix2 = matrix/someValue2;

	ASSERT_DOUBLE_EQ( someValue/someValue2, (matrix2)(0,0) );
}

TEST( Matrix3x3Double, DISABLED_SubtractDivide )
{
	const double someValue = 1.2423;
	const double someValue2 = 4.2332;
	Matrix<double,3> matrix( someValue, someValue, someValue );

	Matrix<double,3> matrix2;

//	TODO: I want to do the following:
//	matrix2 = (matrix-someValue2)/someValue2;

	ASSERT_DOUBLE_EQ( (someValue-someValue2)/someValue2, (matrix2)(0,0) );
}

TEST( Matrix3x3Double, Det )
{
	const double one = 1.;
	const double someValue2 = 4.2332;
	Matrix<double,3> matrix( one, one, one );
	matrix(2,0) = someValue2;
	matrix(0,2) = someValue2;

	ASSERT_DOUBLE_EQ( 1.-someValue2*someValue2, matrix.det() );
}
