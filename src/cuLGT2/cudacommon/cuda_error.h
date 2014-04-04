/**
 * cuda_error.h
 *
 *  Created on: Feb 26, 2014
 *      Author: vogt
 */

#ifndef CUDA_ERROR_H_
#define CUDA_ERROR_H_

#include <stdio.h>

#ifdef NO_CUDA_SAFE_CALL
#define CUDA_SAFE_CALL( call, msg ) call;
#else
#define CUDA_SAFE_CALL( call, msg )\
{\
	cudaError_t ret = call;\
	switch(ret)\
	{\
		case cudaSuccess:\
			 break;\
		default :\
			printf(" ERROR at line %i in %s\n %s: %s\n",\
			__LINE__, __FILE__, msg, cudaGetErrorString(ret));\
			exit(-1);\
			break;\
	}\
}

#define CUDA_LAST_ERROR( msg )\
	{cudaDeviceSynchronize();\
	cudaError_t error = cudaGetLastError();\
	if(error!=cudaSuccess) {\
		fprintf(stderr,"ERROR: %s: %s\n", msg, cudaGetErrorString(error) );\
		exit(-1);\
	}}
#endif

#endif /* CUDA_ERROR_H_ */
