/**
 * Matrix.h
 *
 *  Created on: Jul 4, 2014
 *      Author: vogt
 */

#ifndef MATRIX_H_
#define MATRIX_H_
#include <cmath>
#include <ostream>
#include "cudacommon/cuda_host_device.h"
#include "Vector.h"

namespace culgt
{

template<typename T, int Size> class Matrix
{

};

template<typename T> class Matrix<T,3>
{
public:

	CUDA_HOST_DEVICE Matrix()
	{
		for( int i = 0; i < 9; i++ )
		{
			mat[i] = 0.;
		}
	}

	CUDA_HOST_DEVICE Matrix( T d )
	{
		mat[0] = d;
		mat[1] = 0.;
		mat[2] = 0.;
		mat[3] = 0.;
		mat[4] = d;
		mat[5] = 0.;
		mat[6] = 0.;
		mat[7] = 0.;
		mat[8] = d;
	}

	CUDA_HOST_DEVICE Matrix( T diag0, T diag1, T diag2 )
	{
		for( int i = 0; i < 9; i++ )
		{
			mat[i] = 0.;
		}
		this->operator()( 0, 0 ) = diag0;
		this->operator()( 1, 1 ) = diag1;
		this->operator()( 2, 2 ) = diag2;
	}

	CUDA_HOST_DEVICE Matrix( T* vals )
	{
		for( int i = 0; i < 9; i++ )
		{
			mat[i] = vals[i];
		}
	}

	CUDA_HOST_DEVICE T& operator()( int row, int col )
	{
		return mat[row*3+col];
	}

	CUDA_HOST_DEVICE Vector<T,3> col( int i )
	{
		Vector<T,3> vec;
		for( int j = 0; j < 3; j++ )
		{
			vec[j] = this->operator()(j,i);
		}
		return vec;
	}

	CUDA_HOST_DEVICE void setCol( int i, Vector<T,3> vec )
	{
		for( int j = 0; j < 3; j++ )
		{
			this->operator()(j,i) = vec(j);
		}
	}

	CUDA_HOST_DEVICE Vector<T,3> row( int i )
	{
		Vector<T,3> vec;
		for( int j = 0; j < 3; j++ )
		{
			vec[j] = this->operator()(i,j);
		}
		return vec;
	}

	CUDA_HOST_DEVICE Matrix<T,3>& hermitian()
	{
		for( int i = 0; i < 3; i++ )
		{
			this->operator()(i,i) = conj( this->operator()(i,i) );
			for( int j = 0; j < i; j++ )
			{
				T temp = conj( this->operator()(i,j) );
				this->operator()(i,j) = conj(this->operator()(j,i));
				this->operator()(j,i) = temp;
			}
		}
		return *this;
	}

	CUDA_HOST_DEVICE Matrix<T,3>& operator*=( Matrix<T,3> rhs )
	{
		Matrix<T,3> temp = *this;
//		for( int i = 0; i < 3; i++ )
//			for( int j = 0; j < 3; j++ )
//			{
//				temp(i,j) = this->operator()(i,j);
//			}

		for( int i = 0; i < 3; i++ )
			for( int j = 0; j < 3; j++ )
			{
				this->operator()(i,j) = 0.;
				for( int k = 0; k < 3; k++ )
				{
					this->operator()(i,j) += temp(i,k)*rhs(k,j);
				}
			}
		return *this;
	}

	CUDA_HOST_DEVICE Matrix<T,3>& operator*=( T rhs )
	{
		for( int i = 0; i < 3; i++ )
			for( int j = 0; j < 3; j++ )
			{
				this->operator()(i,j) *= rhs;
			}
		return *this;
	}

	CUDA_HOST_DEVICE Matrix<T,3>& operator/=( T rhs )
	{
		return this->operator*=( 1./rhs );
	}

	CUDA_HOST_DEVICE T trace()
	{
		return this->operator()( 0, 0 )+this->operator()( 1, 1 )+this->operator()( 2, 2 );
	}

	CUDA_HOST_DEVICE T det()
	{
		return ( this->operator()( 0, 0 )*this->operator()( 1, 1 )*this->operator()( 2, 2 )
			+ this->operator()( 0, 1 )*this->operator()( 1, 2 )*this->operator()( 2, 0 )
			+ this->operator()( 0, 2 )*this->operator()( 1, 0 )*this->operator()( 2, 1 )
			- this->operator()( 0, 2 )*this->operator()( 1, 1 )*this->operator()( 2, 0 )
			- this->operator()( 0, 0 )*this->operator()( 1, 2 )*this->operator()( 2, 1 )
			- this->operator()( 0, 1 )*this->operator()( 1, 0 )*this->operator()( 2, 2 ) );
	}

private:
	T mat[9];
};

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator+( Matrix<T,Size>& lhs, T rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		result(i,i) += rhs;
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator-( Matrix<T,Size>& lhs, T rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		result(i,i) -= rhs;
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator*( Matrix<T,Size>& lhs, T rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		for( int j = 0; j < Size; j++ )
		{
			result(i,j) *= rhs;
		}
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator*( Matrix<T,Size>& lhs, Matrix<T,Size>& rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	result *= rhs;
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator/( Matrix<T,Size>& lhs, T rhs )
{
	return lhs*(1./rhs);
}


template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator+( Matrix<T,Size>& lhs, Matrix<T,Size>& rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		for( int j = 0; j < Size; j++ )
			result(i,j) += rhs(i,j);
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Matrix<T,Size> operator-( Matrix<T,Size>& lhs, Matrix<T,Size>& rhs )
{
	Matrix<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		for( int j = 0; j < Size; j++ )
			result(i,j) -= rhs(i,j);
	}
	return result;
}

template<typename T, int Size> inline std::ostream& operator<<(std::ostream& out, Matrix<T,Size>& matrix)
{
	for( int i = 0; i < Size; i++ )
	{
		for( int j = 0; j < Size; j++ )
		{
			out << matrix(i,j) << " \t";
		}
		out << "\n";
	}
	return out;
}

}

#endif /* MATRIX_H_ */
