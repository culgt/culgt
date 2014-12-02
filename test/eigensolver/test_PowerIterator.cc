#include "eigensolver/PowerIterator.h"
#include "gmock/gmock.h"

TEST( TwoXTwoSymmetric, SolvesForDiagonalMatrix )
{
	float matrix[4] = {.5, 0.,0.,1.};
	PowerIterator<float,2> eigen( matrix );

	ASSERT_FLOAT_EQ( 1., eigen.getEigenvalue() );
}

TEST( TwoXTwoSymmetric, SolvesForSimpleMatrix )
{
	float matrix[4] = {1., 2.,3.,4.};
	PowerIterator<float,2> eigen( matrix );

	ASSERT_FLOAT_EQ( 5.372281323269014, eigen.getEigenvalue() );
}

TEST( TwoXTwoSymmetric, EigenvectorForSimpleMatrix )
{
	float matrix[4] = {1., 2.,3.,4.};
	PowerIterator<float,2> eigen( matrix );

	float norm = eigen.getEigenvector()[0]*eigen.getEigenvector()[0]+eigen.getEigenvector()[1]*eigen.getEigenvector()[1];

	float biggerValueInEigenvector = (eigen.getEigenvector()[0]>eigen.getEigenvector()[1])?(eigen.getEigenvector()[0]):(eigen.getEigenvector()[1]);


	ASSERT_FLOAT_EQ( 0.909376709132124, biggerValueInEigenvector/sqrt(norm) );
}

TEST( FourXFourSymmetric, Eigenvalue )
{
	const int SIZE = 4;
	float matrix[SIZE*SIZE] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	PowerIterator<float,SIZE> eigen( matrix, 1e-5 );

	ASSERT_FLOAT_EQ( 36.209372712298553, eigen.getEigenvalue() );
}
