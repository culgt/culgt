/**
 * cuda_gtest_plugin.h
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef CUDA_GTEST_PLUGIN_H_
#define CUDA_GTEST_PLUGIN_H_

#include "gmock/gmock.h"

template <typename T> class type_name {
public:
	typedef T TYPE;
    static const char *name;
};

#define DECLARE_TYPE_NAME(x) template<> const char *type_name<x>::name = #x;
//#define GET_TYPE_NAME(x) (type_name<typeof(x)>::name)
//#define GET_TYPE(x) (typename type_name<typeof(x)>::TYPE)

DECLARE_TYPE_NAME( float );
DECLARE_TYPE_NAME( int );

#define CUDA_ASSERT_FLOAT_EQ(expected,actual)\
	out[0]=expected;\
	out[1]=actual;

#define CUDA_TEST_CLASS_NAME_(test_case_name, test_name)\
  kernel_test_case_name##_##test_name##_Test

#define CUDA_FLOAT_TEST( test_case_name, test_name )\
	__global__ void CUDA_TEST_CLASS_NAME_(test_case_name, test_name)(float*);\
	TEST( test_case_name, test_name )\
	{\
	float* dOut;\
	cudaMalloc( (void**)(&dOut), 2*sizeof( float ) );\
	float out[2] = {0,0};\
	cudaMemcpy( dOut, &out, 2*sizeof(float), cudaMemcpyHostToDevice );\
	CUDA_TEST_CLASS_NAME_(test_case_name, test_name)<<<1,1>>>(dOut);\
	cudaMemcpy( &out, dOut, 2*sizeof(float), cudaMemcpyDeviceToHost );\
	ASSERT_FLOAT_EQ( out[0], out[1] );\
	};\
	__global__ void CUDA_TEST_CLASS_NAME_(test_case_name, test_name)(float* out)


struct TestTransporter
{
	float tfloat[2];
	int tint[2];

	bool evaluateInt;
	bool evaluateFloat;

	TestTransporter()
	{
		evaluateFloat = false;
		tfloat[0] = 0.;
		tfloat[1] = 0.;

		evaluateInt = false;
		tint[0] = 0;
		tint[1] = 0;
	};
};

	template<typename T> static __host__ __device__ void setTestTransporterValue( TestTransporter* transporter, T expected, T actual );

	template<> static __host__ __device__ void setTestTransporterValue( TestTransporter* transporter, float expected, float actual )
	{
		transporter->tfloat[0] = expected;
		transporter->tfloat[1] = actual;
		transporter->evaluateFloat = true;
	}

	template<> static __host__ __device__ void setTestTransporterValue( TestTransporter* transporter, int expected, int actual )
	{
		transporter->tint[0] = expected;
		transporter->tint[1] = actual;
		transporter->evaluateInt = true;
	}

#define CUDA_TEST( test_case_name, test_name )\
	__global__ void CUDA_TEST_CLASS_NAME_(test_case_name, test_name)(TestTransporter*);\
	TEST( test_case_name, test_name )\
	{\
	TestTransporter* dTestTransporter;\
	cudaMalloc( (void**)(&dTestTransporter), sizeof( TestTransporter ) );\
	TestTransporter testTransporter;\
	cudaMemcpy( dTestTransporter, &testTransporter, sizeof(TestTransporter), cudaMemcpyHostToDevice );\
	CUDA_TEST_CLASS_NAME_(test_case_name, test_name)<<<1,1>>>(dTestTransporter);\
	cudaMemcpy( &testTransporter, dTestTransporter, sizeof(TestTransporter), cudaMemcpyDeviceToHost );\
	if( testTransporter.evaluateFloat ) { ASSERT_FLOAT_EQ( testTransporter.tfloat[0], testTransporter.tfloat[1] ); }\
	if( testTransporter.evaluateInt ) { ASSERT_EQ( testTransporter.tint[0], testTransporter.tint[1] ); } \
	};\
	__global__ void CUDA_TEST_CLASS_NAME_(test_case_name, test_name)(TestTransporter* testTransporter)

#define CUDA_ASSERT_EQ(expected,actual)\
		setTestTransporterValue( testTransporter, expected, actual );

#endif /* CUDA_GTEST_PLUGIN_H_ */
