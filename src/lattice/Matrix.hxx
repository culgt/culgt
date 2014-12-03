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
 *
 * Storage class for a NxN matrix. Allocates a N*N linear array with elements of type "T".
 * In contrast to the other storage class "Link", this class actually allocates memory.
 * It has to provide the same functions.
 *
 * TODO:
 *  - see "Link"
 *  - maybe imply symmetric matrix by name.
 */

#ifndef MATRIX_HXX_
#define MATRIX_HXX_

#include "cuda/cuda_host_device.h"
#include "datatype/datatypes.h"
#include <iostream>
#include <fstream>

template<class T, int N> class Matrix
{
public:
	CUDA_HOST_DEVICE inline Matrix();
//	CUDA_HOST_DEVICE inline Matrix( const Matrix<T,N> &m );
	CUDA_HOST_DEVICE inline virtual ~Matrix();
	CUDA_HOST_DEVICE inline T get(int i, int j);
	CUDA_HOST_DEVICE inline void set(int i, int j, T c);
	CUDA_HOST_DEVICE inline T trace();
	CUDA_HOST_DEVICE inline void identity();
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator+=( Matrix<T,N> );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator-=( Matrix<T,N> );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator-=( T );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator/=( T );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator*=( Matrix<T,N> );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator*=( T );
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator=( T );
	CUDA_HOST_DEVICE inline T& operator[]( int i );
	CUDA_HOST_DEVICE inline void print();
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
 * Set to identity.
 */
template<class T, int N> void Matrix<T,N>::identity()
{
	T c = 1; // All types need to support this, especially the complex type!
	for( int i = 0; i < N; i++ )
	{
		set(i,i,c);
	}
}

/**
 * Add and assign...
 */
template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator+=( Matrix<T,N> a )
{
	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			mat[i*N+j] += a.mat[i*N+j];
		}
	}
	return *this;
}

/**
 * Subtract and assign...
 */
template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator-=( Matrix<T,N> a )
{
	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			mat[i*N+j] -= a.mat[i*N+j];
		}
	}
	return *this;
}

template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator-=( T a )
{
	for(int i = 0; i < N; i++ )
	{
		mat[i*N+i] -= a;
	}
	return *this;
}

template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator/=( T a )
{
	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			mat[i*N+j] /= a;
		}
	}
	return *this;
}

/**
 * Multiply and assign...
 * TODO do it better
 */
template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator*=( Matrix<T,N> a )
{
	Matrix<T,N> temp( *this );

	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			mat[i*N+j] = (T)0.;
			for( int k = 0; k < N; k++ )
			{
				mat[i*N+j] += temp.mat[i*N+k] * a.mat[k*N+j];
			}
		}
	}
	return *this;
}

template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator*=( T c )
{
	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			mat[i*N+j] *= c;
		}
	}
	return *this;
}


template<class T, int N> Matrix<T,N>& Matrix<T,N>::operator=( T c )
{
	T zero( 0 );
	for(int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			if( i == j ) set( i, j, c );
			else set( i, j, zero );
		}
	}
	return *this;
}

template<class T, int N> T& Matrix<T,N>::operator[]( int i )
{
	return mat[i];
}

template<class T, int N> void Matrix<T,N>::print()
{
	printf( "[" );
	for( int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			printf( "%f+i*%f", get(i,j).x, get(i,j).y );
			if( j < N-1 ) printf( "\t" );
		}
		if( i < N-1 ) printf( "\n" );
	}
	printf( "]\n" );
}


#endif /* MATRIX_HXX_ */
