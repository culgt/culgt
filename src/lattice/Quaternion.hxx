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
 * Quaternion representation for a SU2 (or proportional -> does not need to be normalized) element.
 *
 * TODO:
 *  - implement all reasonable operations: *,+,...
 *  - This may be a storage class of the not yet implemented "SU2" frontend class. When tackling the question of where to
 *    place the mathematical operations (in frontend or storage classes) consider the comment of the trace() function.
 *  - check if 4-dimension array is compiled different compared to 4 scalar variables.
 *
 * For an introduction to the relation of SU2 and quaternion see for example Gattringer & Lang (2010), p. 81.
 * The first row of a SU2-matrix is given by ( x[0]+i*x[3]    x[2]+i*x[1] ).
 */


#ifndef QUATERNION_HXX_
#define QUATERNION_HXX_


#include "cuda/cuda_host_device.h"
#include "Complex.hxx"
#include "datatype/lattice_typedefs.h"
#include <assert.h>

template<class T> class Quaternion
{
public:
	CUDA_HOST_DEVICE inline Quaternion();
	CUDA_HOST_DEVICE inline virtual ~Quaternion();
	CUDA_HOST_DEVICE inline T& operator[]( int i );
	CUDA_HOST_DEVICE inline Complex<T> get( int i, int j );
	CUDA_HOST_DEVICE inline void set( int i, int j, Complex<T> c );
	CUDA_HOST_DEVICE inline Quaternion<T>& hermitian();
	CUDA_HOST_DEVICE inline T reDet();
	CUDA_HOST_DEVICE inline Complex<T> det();
	CUDA_HOST_DEVICE inline T reTrace();
	CUDA_HOST_DEVICE inline Complex<T> trace();
	CUDA_HOST_DEVICE inline void projectSU2();
	CUDA_HOST_DEVICE inline void zero();
	CUDA_HOST_DEVICE inline void identity();

	CUDA_HOST_DEVICE inline Quaternion<T>& operator*=( Quaternion<T> );
	CUDA_HOST_DEVICE inline Quaternion<T>& operator*=( T );
	CUDA_HOST_DEVICE inline Quaternion<T>& operator+=( Quaternion<T> );
	CUDA_HOST_DEVICE inline Quaternion<T>& operator-=( Quaternion<T> );
	CUDA_HOST_DEVICE inline Quaternion<T>& operator=( Quaternion<T> );


	CUDA_HOST_DEVICE inline void leftMult( Quaternion<T> q );
	CUDA_HOST_DEVICE inline void rightMult( Quaternion<T> q );
private:
	T quat[4];
};

template<class T> Quaternion<T>::Quaternion()
{
}

template<class T> Quaternion<T>::~Quaternion()
{
}

/**
 * Return element of quaternion representation.
 * @parameter index i=0,1,2,3
 * @return element
 */
template<class T> T& Quaternion<T>::operator[]( int i )
{
	return quat[i];
}


/**
 * Returns the element (i,j) if the quaternion is considered as a storage type for SU2.
 * For example get(0,0) returns the complex value x[0]+i*x[3], get(1,1) returns x[0]-i*x[3].
 *
 * TODO element access by this pattern through class "SU2" is not a good idea. -> to many calls to get() for a trace() for example...
 * specialization of SU2 for quaternion type is needed here!!!
 */
template<class T> Complex<T> Quaternion<T>::get( int i, int j )
{
	switch( i )
	{
	case 0:
		switch( j )
		{
		case 0:
			return Complex<T>( quat[0], quat[3] );
		case 1:
			return Complex<T>( quat[2], quat[1] );
		default:
			break;
		}
		break;
	case 1:
		switch( j )
		{
		case 0:
			return Complex<T>( -quat[2], quat[1] );
		case 1:
			return Complex<T>( quat[0], -quat[3] );
		default:
			break;
		}
		break;
	default:
		break;
	}
	assert(false);
	return Complex<T>(0,0);

}

/**
 * TODO see get(); here it is is even more obvious: if we set a SU2 matrix to some value,
 * we do the set in Quaternion twice if we are not aware of the underlying storage class
 */
template<class T> void Quaternion<T>::set( int i, int j, Complex<T> c )
{
	switch( i )
	{
	case 0:
		switch( j )
		{
		case 0:
			quat[0] = c.x;
			quat[3] = c.y;
			break;
		case 1:
			quat[2] = c.x;
			quat[1] = c.y;
			break;
		}
		break;
	case 1:
		switch( j )
		{
		case 0:
			quat[2] = -c.x;
			quat[1] = c.y;
			break;
		case 1:
			quat[0] = c.x;
			quat[3] = -c.y;
			break;
		}
		break;
	}
}

/**
 * Computes the hermitian conjugate.
 */
template<class T> Quaternion<T>& Quaternion<T>::hermitian()
{
	quat[1] = -quat[1];
	quat[2] = -quat[2];
	quat[3] = -quat[3];
	return *this;
}

/**
 * Computes the real part of the determinant
 */
template<class T> T Quaternion<T>::reDet()
{
	return quat[0]*quat[0]+quat[3]*quat[3]+quat[2]*quat[2]+quat[1]*quat[1];
}

/**
 * Computes the det, result as Complex<T> for compatiblity reason (Im det A = 0 always)
 */
template<class T> Complex<T> Quaternion<T>::det()
{
	Complex<T> c( quat[0]*quat[0]+quat[3]*quat[3]+quat[2]*quat[2]+quat[1]*quat[1], 0 );
	return c;
}


/**
 * Computes the trace.
 */
template<class T> T Quaternion<T>::reTrace()
{
	return quat[0]+quat[0];
}

/**
 * Computes the trace.
 */
template<class T> Complex<T> Quaternion<T>::trace()
{
	Complex<T> c( quat[0]+quat[0], 0 );
	return c;
}

/**
 * Projects to SU2
 */
template<class T> void Quaternion<T>::projectSU2()
{
	T det = rsqrt( reDet() );
//	T det = 1./sqrt( reDet() );
	quat[0] *= det;
	quat[1] *= det;
	quat[2] *= det;
	quat[3] *= det;
}

template<class T> void Quaternion<T>::zero()
{
	quat[0] = 0;
	quat[1] = 0;
	quat[2] = 0;
	quat[3] = 0;
}

template<class T> void Quaternion<T>::identity()
{
	quat[0] = 1;
	quat[1] = 0;
	quat[2] = 0;
	quat[3] = 0;
}

template<class T> Quaternion<T>& Quaternion<T>::operator*=( Quaternion<T> q )
{
	T a0 = quat[0]*q[0]-quat[3]*q[3]-quat[1]*q[1]-quat[2]*q[2];
	T a3 = quat[0]*q[3]+quat[3]*q[0]+quat[2]*q[1]-quat[1]*q[2];
	T a2 = quat[0]*q[2]+quat[2]*q[0]-quat[3]*q[1]+quat[1]*q[3];
	T a1 = quat[0]*q[1]+quat[1]*q[0]-quat[2]*q[3]+quat[3]*q[2];
	quat[0] = a0;
	quat[3] = a3;
	quat[2] = a2;
	quat[1] = a1;
	return *this;
}

template<class T> Quaternion<T>& Quaternion<T>::operator*=( T a )
{
	quat[0] *= a;
	quat[3] *= a;
	quat[2] *= a;
	quat[1] *= a;
	return *this;
}

template<class T> Quaternion<T>& Quaternion<T>::operator+=( Quaternion<T> q )
{
	quat[0] += q[0];
	quat[1] += q[1];
	quat[2] += q[2];
	quat[3] += q[3];
	return *this;
}

template<class T> Quaternion<T>& Quaternion<T>::operator-=( Quaternion<T> q )
{
	quat[0] -= q[0];
	quat[1] -= q[1];
	quat[2] -= q[2];
	quat[3] -= q[3];
	return *this;
}

template<class T> Quaternion<T>& Quaternion<T>::operator=( Quaternion<T> q )
{
	quat[0] = q[0];
	quat[1] = q[1];
	quat[2] = q[2];
	quat[3] = q[3];
	return *this;
}

template<class T> void Quaternion<T>::leftMult( Quaternion<T> q )
{
	T a0 = q[0]*quat[0]-q[3]*quat[3]-q[1]*quat[1]-q[2]*quat[2];
	T a3 = q[0]*quat[3]+q[3]*quat[0]+q[2]*quat[1]-q[1]*quat[2];
	T a2 = q[0]*quat[2]+q[2]*quat[0]-q[3]*quat[1]+q[1]*quat[3];
	T a1 = q[0]*quat[1]+q[1]*quat[0]-q[2]*quat[3]+q[3]*quat[2];
	quat[0] = a0;
	quat[3] = a3;
	quat[2] = a2;
	quat[1] = a1;
}



template<class T> void Quaternion<T>::rightMult( Quaternion<T> q )
{
	*this*=q;
}

#endif /* QUATERNION_HXX_ */
