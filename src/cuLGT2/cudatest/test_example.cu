/**
 * test_example.cu
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef TEST_EXAMPLE_CU_
#define TEST_EXAMPLE_CU_

#include "cuda_test_compare.h"
#include "gmock/gmock.h"

struct CudaTestInput
{
	int a;
	int b;
};

struct MyCudaTestOutput
{
	int c;

	bool isEqual( const MyCudaTestOutput& a ) const
	{
		if( c == a.c )
		{
			return true;
		}
		return false;
	}
};


struct test_example_add
{
  __host__ __device__ void operator()( const CudaTestInput& in, MyCudaTestOutput& out ) const
  {
	  out.c = in.a + in.b;
  }
};

TEST( ACudaTest, TestAdd )
{
	CudaTestInput in;
	in.a = 1;
	in.b = 2;

	MyCudaTestOutput out;

	cudaRunAndCompare( test_example_add(), in, out );
}

#endif /* TEST_EXAMPLE_CU_ */
