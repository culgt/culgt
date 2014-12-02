/**
 * test_example.cu
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef TEST_EXAMPLE_CU_
#define TEST_EXAMPLE_CU_

#include "cudatest/cuda_gtest_plugin.h"
#include "gmock/gmock.h"

using namespace std;

CUDA_TEST( ATest, WithFloat )
{
	float expect = 1.23;
	float someVariable = 1.23;

	ASSERT_EQ(expect,someVariable);
}

TEST( ATest, WithFloatOnlyOnHost )
{
	float expect = 1.23;
	float someVariable = 1.23;

	ASSERT_EQ(expect,someVariable);
}

CUDA_TEST( ATest, With2IntThatShouldFailOnFirstAssert )
{
	int expect = 2;
	int someVariable = 3;
	int expect2 = 2;
	int someVariable2 = 2;

	ASSERT_EQ(expect,someVariable);
	ASSERT_EQ(expect2,someVariable2);

}

CUDA_TEST( ATest, With2FloatThatShouldFailOnFirstAssert )
{
	float expect = 2.;
	float someVariable = 3.;
	float expect2 = 2.;
	float someVariable2 = 2.;

	ASSERT_FLOAT_EQ(expect,someVariable);
	ASSERT_FLOAT_EQ(expect2,someVariable2);
}

CUDA_TEST( ATest, ThatShouldFailOnlyOnGPU )
{
	float expect = 1.23;
#ifdef __CUDA_ARCH__
	float someVariable = 1.0;
#else
	float someVariable = 1.23;
#endif
	ASSERT_EQ(expect,someVariable);
}

CUDA_TEST( ATest, WithBool )
{
	bool expect = true;
	bool someVariable = true;

	ASSERT_EQ(expect,someVariable);
}

TEST( ATest, TestThatWeHaveAssertDefinedTheWayWeWantInANonCudaTest )
{
	int x = 1;
	ASSERT_EQ( 1, x );
}


#endif /* TEST_EXAMPLE_CU_ */
