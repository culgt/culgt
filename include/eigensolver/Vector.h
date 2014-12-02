/**
 *
 *  Created on: Jul 4, 2014
 *      Author: vogt
 */

#ifndef VECTOR_H_
#define VECTOR_H_
#include <cmath>
#include <ostream>
#include "cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename T, int Size> class Vector
{

};

//template<typename T> class Vector<Complex<T>,3>
//{
//public:
//
//	CUDA_HOST_DEVICE Vector()
//	{
//		for( int i = 0; i < 3; i++ )
//		{
//			vec[i] = 0.;
//		}
//	}
//
//	CUDA_HOST_DEVICE Complex<T>& operator[](  int i )
//	{
//		return vec[i];
//	}
//	CUDA_HOST_DEVICE Complex<T>& operator()(  int i )
//	{
//		return vec[i];
//	}
//
//	CUDA_HOST_DEVICE T norm()
//	{
//		return sqrt( vec[0].abs_squared() +vec[1].abs_squared()+vec[2].abs_squared());
//	}
//
//
//private:
//	Complex<T> vec[3];
//};



template<typename T> class Vector<T,3>
{
public:

	CUDA_HOST_DEVICE Vector()
	{
		for( int i = 0; i < 3; i++ )
		{
			vec[i] = 0.;
		}
	}

	CUDA_HOST_DEVICE T& operator[](  int i )
	{
		return vec[i];
	}
	CUDA_HOST_DEVICE T& operator()(  int i )
	{
		return vec[i];
	}

	CUDA_HOST_DEVICE T norm()
	{
		return sqrt( norm_squared( vec[0] ) +norm_squared( vec[1])+norm_squared( vec[2]));
	}

	CUDA_HOST_DEVICE Vector<T,3>& operator-=( Vector<T,3> rhs )
	{
		for( int i = 0; i < 3; i++ )
		{
				this->operator()(i) -= rhs(i);
		}
		return *this;
	}


private:
	T vec[3];
};


template<typename T, int Size> CUDA_HOST_DEVICE  Vector<T,Size> operator*( Vector<T,Size>& lhs, T rhs )
{
	Vector<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		result(i) *= rhs;
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE  Vector<T,Size> operator/( Vector<T,Size>& lhs, T rhs )
{
	return lhs*(1./rhs);
}

template<typename T, int Size> CUDA_HOST_DEVICE  Vector<T,Size> operator-( Vector<T,Size>& lhs, Vector<T,Size>& rhs )
{
	Vector<T,Size> result;
	result = lhs;
	for( int i = 0; i < Size; i++ )
	{
		result(i) -= rhs(i);
	}
	return result;
}

template<typename T, int Size> CUDA_HOST_DEVICE inline std::ostream& operator<<(std::ostream& out, Vector<T,Size>& vec)
{
	for( int i = 0; i < Size; i++ )
	{
			out << vec(i) << "\n";
	}
	return out;
}

}

#endif
