/*
 * PowerIterator.h
 *
 *  Created on: 04.04.2014
 *      Author: vogt
 */

#ifndef POWERITERATOR_H_
#define POWERITERATOR_H_

#include <cmath>
#include "cudacommon/cuda_host_device.h"
using namespace std;

template<typename T, int Size> class PowerIterator
{
public:
	CUDA_HOST_DEVICE inline PowerIterator( T* mat, T precision = 1e-6, int maxiter = 15 ) : mat(mat), precision(precision), maxiter(maxiter)
	{
		compute();
	}

	CUDA_HOST_DEVICE inline T getEigenvalue()
	{
		return eigenvalue;
	}

	CUDA_HOST_DEVICE inline T* getEigenvector()
	{
		return eigenvector;
	}
private:
	CUDA_HOST_DEVICE inline void compute()
	{
		T dd = 1.0;
		for( int i = 0; i < Size; i++ )
		{
			eigenvector[i] = .5;
		}
	   	T lambda = 10.0;
	   	T y[Size];

	   	int iter = 0;
	   	for( iter = 0; iter < maxiter; iter++ )
//	   	while( dd > precision )
	   	{
			mult( y, mat, eigenvector ); //y = A*x;
			dd = fabs(norm(eigenvector) - lambda);
			lambda = norm(eigenvector);
			div(eigenvector, y, lambda);//x = y / lambda;
			if( dd < precision ) break;
	   }
	   eigenvalue = lambda;
	   div( eigenvector, eigenvector, norm( eigenvector ) );
	}

	CUDA_HOST_DEVICE inline void mult( T* dest, T* mat, T* vec )
	{
		for( int i = 0; i < Size; i++ )
		{
			dest[i] = 0;
			for( int j = 0; j < Size; j++ )
			{
				dest[i] += mat[i*Size+j]*vec[j];
			}
		}
	}

	CUDA_HOST_DEVICE inline void div( T* dest, T* vec, T scalar )
	{
		for( int i = 0; i < Size; i++ )
		{
			dest[i] = vec[i] / scalar;
		}
	}

	CUDA_HOST_DEVICE inline T norm( T* vec )
	{
		T result = 0.;
		for( int i = 0; i < Size; i++ )
		{
			result += vec[i]*vec[i];
		}
		return sqrt(result);
	}

	CUDA_HOST_DEVICE inline T fabs( T val )
	{
		return (val < 0 )?(-val):(val);
	}

	T* mat;
	T eigenvalue;
	T eigenvector[Size];
	T precision;
	int maxiter;
};


#endif /* POWERITERATOR_H_ */
