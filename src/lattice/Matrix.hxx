/**
 * Storage class for a NxN matrix. Allocates a N*N linear array with elements of type "T".
 * In contrast to the other storage class "Link", this class actually allocates memory.
 * It has to provide the same functions.
 *
 * TODO:
 *  - see "Link"
 */

#ifndef MATRIX_HXX_
#define MATRIX_HXX_

#include "../util/datatype/datatypes.h"
#include <iostream>

template<class T, int N> class Matrix
{
public:
	CUDA_HOST_DEVICE inline Matrix();
	CUDA_HOST_DEVICE inline virtual ~Matrix();
	CUDA_HOST_DEVICE inline T get(int i, int j);
	CUDA_HOST_DEVICE inline void set(int i, int j, T c);
	CUDA_HOST_DEVICE inline T trace();
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator+=( Matrix<T,N> );
	T mat[N*N]; // array keeping the matrix TODO make this private?
private:
};

template<class T, int N> Matrix<T,N>::Matrix()
{
}

template<class T, int N> Matrix<T,N>::~Matrix()
{
}

/**
 * Returns the matrix element (i,j).
 * @parameter row index i
 * @parameter col index j
 * @return element (i,j)
 */
template<class T, int N> T Matrix<T,N>::get( int i, int j )
{
	return mat[i*N+j];
}

/**
 * Sets the matrix element (i,j).
 * @parameter row index i
 * @parameter col index j
 * @parameter element to set
 */
template<class T, int N> void Matrix<T,N>::set( int i, int j, T c )
{
	mat[i*N+j] = c;
}

/**
 * Trace.
 */
template<class T, int N> T Matrix<T,N>::trace()
{
	T c;
	for( int i = 0; i < N; i++ )
	{
		c += get(i,i);
	}
	return c;
}

/**
 * Add and assign...
 */
template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator+=( Matrix<T,N> a )
{
	for(int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			mat[i*N+j] += a.mat[i*N+j];
		}
	}
	return *this;
}



#endif /* MATRIX_HXX_ */
