#include "gmock/gmock.h"
#include "util/reduction/Reduction.h"
#include "cudacommon/DeviceCommunicator.h"
#include "math/Complex.h"

using namespace testing;
using namespace culgt;
using namespace std;

TEST( AReduction, DoubleReduction512 )
{
	const int SIZE = 512;
	double* data;
	CUDA_SAFE_CALL( cudaMalloc( &data, SIZE*sizeof(double) ), "Field: Allocating Device memory" );

	for( int i = 0; i < SIZE; i++ )
	{
		DeviceCommunicator<double>::setValue( data, i, 1.0 );
	}

	Reduction<double> reduction( SIZE );

	ASSERT_DOUBLE_EQ( 512., reduction.reduceAll( data ) );
}

TEST( AReduction, DoubleReductionDot512 )
{
	const int SIZE = 512;
	double* data;
	CUDA_SAFE_CALL( cudaMalloc( &data, SIZE*sizeof(double) ), "Field: Allocating Device memory" );
	double* data2;
	CUDA_SAFE_CALL( cudaMalloc( &data2, SIZE*sizeof(double) ), "Field: Allocating Device memory" );

	for( int i = 0; i < SIZE; i++ )
	{
		DeviceCommunicator<double>::setValue( data, i, 1.0 );
		DeviceCommunicator<double>::setValue( data2, i, 2.0 );
	}

	Reduction<double> reduction( SIZE );

	ASSERT_DOUBLE_EQ( 1024., reduction.reduceAllDot( data, data2 ) );
}

TEST( AReduction, ComplexDoubleReduction512 )
{
	const int SIZE = 512;
	Complex<double>* data;
	CUDA_SAFE_CALL( cudaMalloc( &data, SIZE*sizeof(Complex<double>) ), "Field: Allocating Device memory" );

	for( int i = 0; i < SIZE; i++ )
	{
		DeviceCommunicator<Complex<double> >::setValue( data, i, Complex<double>(1.0, 0) );
	}

	Reduction<Complex<double> > reduction( SIZE );

	ASSERT_DOUBLE_EQ( 512., reduction.reduceAll( data ).x );
}

TEST( AReduction, ComplexDoubleReductionDotConjugate512 )
{
	const int SIZE = 512;
	Complex<double>* data;
	CUDA_SAFE_CALL( cudaMalloc( &data, SIZE*sizeof(Complex<double>) ), "Field: Allocating Device memory" );
	Complex<double>* data2;
	CUDA_SAFE_CALL( cudaMalloc( &data2, SIZE*sizeof(Complex<double>) ), "Field: Allocating Device memory" );

	for( int i = 0; i < SIZE; i++ )
	{
		DeviceCommunicator<Complex<double> >::setValue( data, i, Complex<double>( 1.0, 0.0) );
		DeviceCommunicator<Complex<double> >::setValue( data2, i, Complex<double>( 0.0, -2.0) );
	}

	Reduction<Complex<double> > reduction( SIZE );

	ASSERT_DOUBLE_EQ( 1024., reduction.reduceAllDotConjugate( data, data2 ).y );
}

TEST( AReduction, DotConjugateIsCompatibleWithDouble512 )
{
	const int SIZE = 512;
	double* data;
	CUDA_SAFE_CALL( cudaMalloc( &data, SIZE*sizeof(double) ), "Field: Allocating Device memory" );
	double* data2;
	CUDA_SAFE_CALL( cudaMalloc( &data2, SIZE*sizeof(double) ), "Field: Allocating Device memory" );

	for( int i = 0; i < SIZE; i++ )
	{
		DeviceCommunicator<double>::setValue( data, i, 1. );
		DeviceCommunicator<double>::setValue( data2, i, 2. );
	}

	Reduction<double> reduction( SIZE );

	ASSERT_DOUBLE_EQ( 1024., reduction.reduceAllDotConjugate( data, data2 ) );
}
