/**
 * cuda_test_outputtypes.h
 *
 *  Created on: Feb 24, 2014
 *      Author: vogt
 */

#ifndef CUDA_TEST_OUTPUTTYPES_H_
#define CUDA_TEST_OUTPUTTYPES_H_


struct CudaTestOutputFloat
{
	float result;

	bool isEqual( CudaTestOutputFloat& a )
	{
		if( a.result == result ) return true;
		else return false;
	}
};



#endif /* CUDA_TEST_OUTPUTTYPES_H_ */
