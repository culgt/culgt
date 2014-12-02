/**
 * test_example.cu
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef TEST_EXAMPLE_GTEST_PLUGIN_HOST_CU_
#define TEST_EXAMPLE_GTEST_PLUGIN_HOST_CU_

#include "cudatest/cuda_gtest_plugin.h"
#include "gmock/gmock.h"

using namespace std;

CUDA_TEST( ATest, WithFloatNotCompiledWithCudaShouldGiveWarningOnCompile )
{
	float expect = 1.23;
	float someVariable = 1.23;

	ASSERT_EQ(expect,someVariable);
}

#endif /* TEST_EXAMPLE_CU_ */
