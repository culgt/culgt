/*
 * CudaError.hxx
 *
 *  Created on: Nov 12, 2012
 *      Author: vogt
 */

#ifndef CUDAERROR_HXX_
#define CUDAERROR_HXX_

#include <iostream>
#include <cuda.h>

class CudaError
{
public:
	static void getLastError( const char *msg );
	static void safeCall( cudaError_t call, const char *msg = "" );
};


void CudaError::getLastError( const char *msg )
{
	cudaError_t error = cudaGetLastError();
	if(error!=cudaSuccess) {
		fprintf(stderr,"ERROR: %s: %s\n", msg, cudaGetErrorString(error) );
		exit(-1);
	}

}


void CudaError::safeCall( cudaError_t call, const char *msg )
{

	cudaError_t ret = call;
	switch(ret)
	{
		case cudaSuccess:
			 break;
		default :
			printf(" ERROR at line : %i\n %s: %s\n",
			__LINE__, msg, cudaGetErrorString(ret));
			exit(-1);
			break;
	}
}

#endif /* CUDAERROR_HXX_ */
