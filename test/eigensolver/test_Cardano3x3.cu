#include "gmock/gmock.h"
#include <iostream>
#include "math/Complex.h"
#include "cudatest/cuda_gtest_plugin.h"
#include "eigensolver/Cardano3x3.h"

using namespace culgt;
using namespace std;



//CUDA_TEST( Cardano3x3Float, ComputesEigenvaluesForFullMatrix )
//{
//	float values[9] = {0.5648000,   0.9337000 ,  0.5146000,
//			   0.0044000 ,  0.7516000 ,  0.0349000,
//			   0.9919000,   0.7123000 ,  0.7573000 };
//	float eigenvals[3] = {    1.438722545459362e+00,
//		    -4.629879539880649e-02,
//		     6.812762499394441e-01};
//
//	Matrix<float,3> matrix( values );
//
//	Cardano3x3<float> eigen( matrix );
//
//	ASSERT_FLOAT_EQ( eigenvals[0], eigen.getEigenvalue( 0 ) );
//	ASSERT_FLOAT_EQ( eigenvals[1], eigen.getEigenvalue( 1 ) );
//	ASSERT_FLOAT_EQ( eigenvals[2], eigen.getEigenvalue( 2 ) );
//}
