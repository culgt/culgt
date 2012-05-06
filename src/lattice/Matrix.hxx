/*
 * Link.hxx
 *
 *  Created on: Apr 17, 2012
 *      Author: vogt
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
//	CUDA_HOST_DEVICE inline complex get(int iSub, int jSub, int i, int j);
	CUDA_HOST_DEVICE inline void set(int i, int j, T c);
//	CUDA_HOST_DEVICE inline void set(int iSub, int jSub, int i, int j, complex c);
	CUDA_HOST_DEVICE inline T trace();
//	CUDA_HOST_DEVICE inline complex det();
	CUDA_HOST_DEVICE inline Matrix<T, N>& operator+=( Matrix<T,N> );
	T mat[N*N];
private:
};

template<class T, int N> Matrix<T,N>::Matrix()
{
}

template<class T, int N> Matrix<T,N>::~Matrix()
{
}


template<class T, int N> T Matrix<T,N>::get( int i, int j )
{
	return mat[i*N+j];
}

//template<class Pattern, class TheSite, int T_Ndim, int T_Nc> complex Link<Pattern, TheSite, T_Ndim, T_Nc>::get( int i, int j )
//{
//	// TODO check if make_cucomplex is overloaded for double precision usage.
//	return complex( data[Pattern::getIndex( site, mu, i, j, 0 )], data[Pattern::getIndex( site, mu, i, j, 1 )] );
//}

template<class T, int N> void Matrix<T,N>::set( int i, int j, T c )
{
	mat[i*N+j] = c;
}

template<class T, int N> T Matrix<T,N>::trace()
{
	T c;
	for( int i = 0; i < N; i++ )
	{
		c += get(i,i);
	}
	return c;
}

//template<class Pattern, class TheSite, int T_Ndim, int T_Nc> complex Link<Pattern, TheSite, T_Ndim, T_Nc>::det()
//{
//	assert( T_Nc == 3 ); // TODO do it properly
//	complex c( 0, 0 );
//	complex temp( 1, 0 );
//	temp *= get(0,0);
//	temp *= get(1,1);
//	temp *= get(2,2);
//	c += temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,1);
//	temp *= get(1,2);
//	temp *= get(2,0);
//	c += temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,2);
//	temp *= get(1,0);
//	temp *= get(2,1);
//	c += temp;
//
//
//
//	temp = complex( 1, 0 );
//	temp *= get(2,0);
//	temp *= get(1,1);
//	temp *= get(0,2);
//	c -= temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(1,0);
//	temp *= get(0,1);
//	temp *= get(2,2);
//	c -= temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,0);
//	temp *= get(2,1);
//	temp *= get(1,2);
//	c -= temp;
//
//	return c;
//}


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
