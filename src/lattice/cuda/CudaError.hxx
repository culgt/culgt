/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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
