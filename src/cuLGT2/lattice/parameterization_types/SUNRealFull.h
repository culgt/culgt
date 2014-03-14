/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SUNREALFULL_H_
#define SUNREALFULL_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"

namespace culgt
{

/*
 *
 */
template<int Nc, typename T> class SUNRealFull
{
public:
	static const lat_dim_t NC = Nc;
	static const lat_group_index_t SIZE = 2*(Nc*Nc);
	typedef T TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
		for( int i = 0; i < Nc; i++ )
		{
			store[2*i*(1+Nc)] = 1.0;
		}
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		REALTYPE result = 0;
		for( int i = 0; i < Nc; i++ )
		{
			result += store[2*i*(1+Nc)];
		}
		return result;
	}
};

template<typename T> class SUNRealFull<3,T>
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = 2*(3*3);
	typedef T TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
		for( int i = 0; i < 3; i++ )
		{
			store[2*i*(1+3)] = 1.0;
		}
	}

	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
	{
		return store[0]*store[8]*store[16] - store[0]*store[10]*store[14] - store[2]*store[6]*store[16] + store[2]*store[10]*store[12] + store[4]*store[6]*store[14] - store[4]*store[8]*store[12] - store[0]*store[9]*store[17] + store[0]*store[11]*store[15] - store[1]*store[8]*store[17] - store[1]*store[9]*store[16] + store[1]*store[10]*store[15] + store[1]*store[11]*store[14] + store[2]*store[7]*store[17] - store[2]*store[11]*store[13] + store[3]*store[6]*store[17] + store[3]*store[7]*store[16] - store[3]*store[10]*store[13] - store[3]*store[11]*store[12] - store[4]*store[7]*store[15] + store[4]*store[9]*store[13] - store[5]*store[6]*store[15] - store[5]*store[7]*store[14] + store[5]*store[8]*store[13] + store[5]*store[9]*store[12];
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return store[0]+store[8]+store[16];
	}

	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		TYPE a[SIZE];
		for( int i = 0; i < SIZE; i++ )
		{
			a[i] = dest[i];
		}

		dest[0] = a[0]*b[0] - a[1]*b[1] + a[2]*b[6] - a[3]*b[7] + a[4]*b[12] - a[5]*b[13];
		dest[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[7] + a[3]*b[6] + a[4]*b[13] + a[5]*b[12];
		dest[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[8] - a[3]*b[9] + a[4]*b[14] - a[5]*b[15];
		dest[3] = a[0]*b[3] + a[1]*b[2] + a[2]*b[9] + a[3]*b[8] + a[4]*b[15] + a[5]*b[14];
		dest[4] = a[0]*b[4] - a[1]*b[5] + a[2]*b[10] - a[3]*b[11] + a[4]*b[16] - a[5]*b[17];
		dest[5] = a[0]*b[5] + a[1]*b[4] + a[2]*b[11] + a[3]*b[10] + a[4]*b[17] + a[5]*b[16];
		dest[6] = a[6]*b[0] - a[7]*b[1] + a[8]*b[6] - a[9]*b[7] + a[10]*b[12] - a[11]*b[13];
		dest[7] = a[6]*b[1] + a[7]*b[0] + a[8]*b[7] + a[9]*b[6] + a[10]*b[13] + a[11]*b[12];
		dest[8] = a[6]*b[2] - a[7]*b[3] + a[8]*b[8] - a[9]*b[9] + a[10]*b[14] - a[11]*b[15];
		dest[9] = a[6]*b[3] + a[7]*b[2] + a[8]*b[9] + a[9]*b[8] + a[10]*b[15] + a[11]*b[14];
		dest[10] = a[6]*b[4] - a[7]*b[5] + a[8]*b[10] - a[9]*b[11] + a[10]*b[16] - a[11]*b[17];
		dest[11] = a[6]*b[5] + a[7]*b[4] + a[8]*b[11] + a[9]*b[10] + a[10]*b[17] + a[11]*b[16];
		dest[12] = a[12]*b[0] - a[13]*b[1] + a[14]*b[6] - a[15]*b[7] + a[16]*b[12] - a[17]*b[13];
		dest[13] = a[12]*b[1] + a[13]*b[0] + a[14]*b[7] + a[15]*b[6] + a[16]*b[13] + a[17]*b[12];
		dest[14] = a[12]*b[2] - a[13]*b[3] + a[14]*b[8] - a[15]*b[9] + a[16]*b[14] - a[17]*b[15];
		dest[15] = a[12]*b[3] + a[13]*b[2] + a[14]*b[9] + a[15]*b[8] + a[16]*b[15] + a[17]*b[14];
		dest[16] = a[12]*b[4] - a[13]*b[5] + a[14]*b[10] - a[15]*b[11] + a[16]*b[16] - a[17]*b[17];
		dest[17] = a[12]*b[5] + a[13]*b[4] + a[14]*b[11] + a[15]*b[10] + a[16]*b[17] + a[17]*b[16];
	}

	static CUDA_HOST_DEVICE void inline flipSign( TYPE& a )
	{
		a = -a;
	}

	static CUDA_HOST_DEVICE void inline swap( TYPE& a, TYPE& b )
	{
		TYPE temp = b;
		b = a;
		a = temp;
	}

	static CUDA_HOST_DEVICE void inline swapAndFlipSign( TYPE& a, TYPE& b )
	{
		TYPE temp = b;
		b = -a;
		a = -temp;
	}

	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
	{
		flipSign( store[1] );
		swap( store[2], store[6] );
		swapAndFlipSign( store[3], store[7] );
		swap( store[4], store[12] );
		swapAndFlipSign( store[5], store[13] );
		flipSign( store[9] );
		swap( store[10], store[14] );
		swapAndFlipSign( store[11], store[15] );
		flipSign( store[17] );
	}
};


} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
