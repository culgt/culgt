/*
 * Quaternion.hxx
 *
 *  Created on: Apr 20, 2012
 *      Author: vogt
 */

#ifndef QUATERNION_HXX_
#define QUATERNION_HXX_

template<class T> class Quaternion
{
public:
	CUDA_HOST_DEVICE inline Quaternion();
	CUDA_HOST_DEVICE inline virtual ~Quaternion();
	CUDA_HOST_DEVICE inline T& operator[]( int i );
	CUDA_HOST_DEVICE inline Complex<T> get( int i, int j );
	CUDA_HOST_DEVICE inline void set( int i, int j, Complex<T> c );
	CUDA_HOST_DEVICE inline Quaternion<T>& hermitian();
	CUDA_HOST_DEVICE inline Complex<T> det();
private:
	T quat[4];
};

template<class T> Quaternion<T>::Quaternion()
{
}

template<class T> Quaternion<T>::~Quaternion()
{
}

template<class T> T& Quaternion<T>::operator[]( int i )
{
	return quat[i];
}


/**
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
		}
	case 1:
		switch( j )
		{
		case 0:
			return Complex<T>( -quat[2], quat[1] );
		case 1:
			return Complex<T>( quat[0], -quat[3] );
		}
	}
	return Complex<T>(0,0);
}

/**
 * TODO see get(); here it is is even more obvious: if we set a SU2 matrix to some value we do the set in Quaternion twice
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
		}
	}
}

template<class T> Quaternion<T>& Quaternion<T>::hermitian()
{
//	Quaternion<T> q;
//	q[0] = quat[0];
//	q[1] = -quat[1];
//	q[2] = -quat[2];
//	q[3] = -quat[3];

	quat[1] = -quat[1];
	quat[2] = -quat[2];
	quat[3] = -quat[3];
	return *this;
}

template<class T> Complex<T> Quaternion<T>::det()
{
	Complex<T> c( quat[0]*quat[0]-quat[3]*quat[3]+quat[2]*quat[2]-quat[1]*quat[1], 0 ); // TODO im
	return c;
}

#endif /* QUATERNION_HXX_ */
