/**
 * test_example.cu
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef TEST_EXAMPLE_CU_
#define TEST_EXAMPLE_CU_

#include "cuda_gtest_plugin.h"
#include "gmock/gmock.h"
#include <iostream>

using namespace std;


CUDA_FLOAT_TEST( ACudaTest, Test )
{
	float expect = 1.23;
	float someVariable = 1.23;

	CUDA_ASSERT_FLOAT_EQ(expect,someVariable);
};



TEST( ATest, Test )
{
	// check if I can find out the type of a variable
	float x = 1.23;

	typename type_name<typeof(x)>::TYPE y = x;

	ASSERT_EQ( x, y );
}


CUDA_TEST( ATest, Test2 )
{
	float expect = 1.23;
	float someVariable = 1.23;

	CUDA_ASSERT_EQ(expect,someVariable);
}

CUDA_TEST( ATest, Test3 )
{
	int expect = 2;
	int someVariable = 2;

	CUDA_ASSERT_EQ(expect,someVariable);
}

CUDA_TEST( ATest, DoesNotCareIfIMessWithTransporter )
{
	int expect = 2;
	int someVariable = 2;

	testTransporter->tfloat[0] = 1.;
	testTransporter->tfloat[1] = 2.;

	CUDA_ASSERT_EQ(expect,someVariable);
}


/*
 *  DO NOT DELETE THIS: IT IS FOR ILLUSTRATION:
 */
//CUDA_TEST( ATest, DoesNotCompileIfIUseNonCudaAssert )
//{
//	int expect = 2;
//	int someVariable = 2;
//	ASSERT_EQ(expect,someVariable);
//}
/*
 * END: DO NOT DELETE THIS
 */




#endif /* TEST_EXAMPLE_CU_ */
