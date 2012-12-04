/*
 * Complex.hxx
 *
 *  Created on: Apr 18, 2012
 *      Author: vogt
 */

#ifndef COMPLEX_HXX_
#define COMPLEX_HXX_

#include "cuda/cuda_host_device.h"

template <class datatype> class Complex
{
public:
	CUDA_HOST_DEVICE inline Complex( const datatype x, const datatype y );
	CUDA_HOST_DEVICE inline Complex( const datatype x );
	CUDA_HOST_DEVICE inline Complex();
	CUDA_HOST_DEVICE inline virtual ~Complex();
	CUDA_HOST_DEVICE inline datatype abs();
	CUDA_HOST_DEVICE inline datatype abs_squared();
	CUDA_HOST_DEVICE inline Complex<datatype> conj();
	CUDA_HOST_DEVICE inline Complex<datatype>& operator+=( const Complex<datatype> a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator+=( const datatype a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator-=( const Complex<datatype> a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator-=( const datatype a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator*=( const Complex<datatype> a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator*=( const datatype a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator/=( const datatype a );
	CUDA_HOST_DEVICE inline Complex<datatype>& operator/=( const Complex<datatype> a );
	datatype x; // we use the naming convention of cuComplex.h
	datatype y;
};



template<class datatype> Complex<datatype>::Complex(const datatype x, const datatype y)
{
	this->x = x;
	this->y = y;
}

template<class datatype> Complex<datatype>::Complex(const datatype x )
{
	this->x = x;
	this->y = 0;
}

template<class datatype> Complex<datatype>::Complex()
{
	this->x = 0;
	this->y = 0;
}


template<class datatype> Complex<datatype>::~Complex()
{
}


template<class datatype> datatype Complex<datatype>::abs()
{
	return sqrt( x*x + y*y );
}

template<class datatype> datatype Complex<datatype>::abs_squared()
{
	return ( x*x + y*y );
}

template<class datatype> Complex<datatype> Complex<datatype>::conj()
{
	return Complex<datatype>(x,-y);
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator+=( const Complex<datatype> a )
{
	x += a.x;
	y += a.y;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator+=( const datatype a )
{
	x += a;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator-=( const Complex<datatype> a )
{
	x -= a.x;
	y -= a.y;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator-=( const datatype a )
{
	x -= a;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator*=( const Complex<datatype> a )
{
	datatype re = x;
	datatype im = y;
	x = re*a.x - im*a.y;
	y = re*a.y + im*a.x;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator*=( const datatype a )
{
	x *= a;
	y *= a;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator/=( const datatype a )
{
	x /= a;
	y /= a;
	return *this;
}

template<class datatype> Complex<datatype>& Complex<datatype>::operator/=( const Complex<datatype> a )
{
	datatype re = x;
	datatype im = y;
	datatype abs = a.x*a.x+a.y*a.y;
	x = (re*a.x + im*a.y)/abs;
	y = (-re*a.y + im*a.x)/abs;

	return *this;
}







template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator+( Complex<datatype> a, Complex<datatype> b )
{
	Complex<datatype> c = a;
	return c+=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator+( Complex<datatype> a, datatype b )
{
	Complex<datatype> c = a;
	return c+=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator+( datatype a, Complex<datatype> b )
{
	Complex<datatype> c(a);
	return c+=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator-( Complex<datatype> a, Complex<datatype> b )
{
	Complex<datatype> c = a;
	return c-=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator-( Complex<datatype> a, datatype b )
{
	Complex<datatype> c = a;
	return c-=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator-( datatype a, Complex<datatype> b )
{
	Complex<datatype> c(a);
	return c-=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator*( Complex<datatype> a, Complex<datatype> b )
{
	Complex<datatype> c = a;
	return c*=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator*( Complex<datatype> a, datatype b )
{
	Complex<datatype> c = a;
	return c*=b;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator*( datatype a, Complex<datatype> b )
{
	Complex<datatype> c = b; // we commuted arguments to use *=(datatype) which is simpler
	return c*=a;
}

template<class datatype> CUDA_HOST_DEVICE static Complex<datatype> operator/( Complex<datatype> a, datatype b )
{
	Complex<datatype> c = a;
	return c/=b;
}

#endif /* COMPLEX_HXX_ */
