/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SUNCOMPLEXFULL_H_
#define SUNCOMPLEXFULL_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"
#include "../../cuLGT1legacy/Complex.hxx"

namespace culgt
{

/*
 *
 */
template<int Nc, typename T> class SUNComplexFull
{
public:
	static const lat_dim_t NC = Nc;
	static const lat_group_index_t SIZE = (Nc*Nc);
	typedef Complex<T> TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = TYPE( 0, 0 );
		}
		for( int i = 0; i < Nc; i++ )
		{
			store[i*(1+Nc)] = TYPE( 1, 0 );
		}
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		REALTYPE result = 0;
		for( int i = 0; i < Nc; i++ )
		{
			result += store[i*(1+Nc)].x;
		}
		return result;
	}

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = TYPE( 0, 0 );
		}
	}
};

//template<typename T> class SUNRealFull<2,T>
//{
//public:
//	static const lat_dim_t NC = 2;
//	static const lat_group_index_t SIZE = 8;
//	typedef T TYPE;
//	typedef T REALTYPE;
//
//	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
//	{
//		for( int i = 0; i < SIZE; i++ )
//		{
//			store[i] = 0;
//		}
//	}
//
//	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
//	{
//		store[0] = 1.;
//		store[1] = 0.;
//		store[2] = 0.;
//		store[3] = 0.;
//		store[4] = 0.;
//		store[5] = 0.;
//		store[6] = 1.;
//		store[7] = 0.;
//	}
//
//	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
//	{
//		return store[0]*store[6] - store[2]*store[4] - store[1]*store[7] + store[3]*store[5];
//	}
//
//	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
//	{
//		return store[0]+store[6];
//	}
//
//	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
//	{
//		TYPE a[SIZE];
//		for( int i = 0; i < SIZE; i++ )
//		{
//			a[i] = dest[i];
//		}
//
//		dest[0] = a[0]*b[0] - a[1]*b[1] + a[2]*b[4] - a[3]*b[5];
//		dest[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[5] + a[3]*b[4];
//		dest[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[6] - a[3]*b[7];
//		dest[3] = a[0]*b[3] + a[1]*b[2] + a[2]*b[7] + a[3]*b[6];
//		dest[4] = a[4]*b[0] - a[5]*b[1] + a[6]*b[4] - a[7]*b[5];
//		dest[5] = a[4]*b[1] + a[5]*b[0] + a[6]*b[5] + a[7]*b[4];
//		dest[6] = a[4]*b[2] - a[5]*b[3] + a[6]*b[6] - a[7]*b[7];
//		dest[7] = a[4]*b[3] + a[5]*b[2] + a[6]*b[7] + a[7]*b[6];
//	}
//
//	static CUDA_HOST_DEVICE void inline flipSign( TYPE& a )
//	{
//		a = -a;
//	}
//
//	static CUDA_HOST_DEVICE void inline swap( TYPE& a, TYPE& b )
//	{
//		TYPE temp = b;
//		b = a;
//		a = temp;
//	}
//
//	static CUDA_HOST_DEVICE void inline swapAndFlipSign( TYPE& a, TYPE& b )
//	{
//		TYPE temp = b;
//		b = -a;
//		a = -temp;
//	}
//
//	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
//	{
//		flipSign( store[1] );
//		swap( store[2], store[4] );
//		swapAndFlipSign( store[3], store[5] );
//		flipSign( store[7] );
//	}
//};

template<typename T> class SUNComplexFull<3,T>
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = (3*3);
	typedef Complex<T> TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = TYPE( 0, 0 );
		}
	}

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = TYPE( 0, 0 );
		}
		for( int i = 0; i < NC; i++ )
		{
			store[i*(1+NC)] = TYPE( 1, 0 );
		}
	}

//	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
//	{
//		return store[0]*store[8]*store[16] - store[0]*store[10]*store[14] - store[2]*store[6]*store[16] + store[2]*store[10]*store[12] + store[4]*store[6]*store[14] - store[4]*store[8]*store[12] - store[0]*store[9]*store[17] + store[0]*store[11]*store[15] - store[1]*store[8]*store[17] - store[1]*store[9]*store[16] + store[1]*store[10]*store[15] + store[1]*store[11]*store[14] + store[2]*store[7]*store[17] - store[2]*store[11]*store[13] + store[3]*store[6]*store[17] + store[3]*store[7]*store[16] - store[3]*store[10]*store[13] - store[3]*store[11]*store[12] - store[4]*store[7]*store[15] + store[4]*store[9]*store[13] - store[5]*store[6]*store[15] - store[5]*store[7]*store[14] + store[5]*store[8]*store[13] + store[5]*store[9]*store[12];
//	}
//
//	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
//	{
//		return store[0]+store[8]+store[16];
//	}
//
//	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
//	{
//		TYPE a[SIZE];
//		for( int i = 0; i < SIZE; i++ )
//		{
//			a[i] = dest[i];
//		}
//
//		dest[0] = a[0]*b[0] - a[1]*b[1] + a[2]*b[6] - a[3]*b[7] + a[4]*b[12] - a[5]*b[13];
//		dest[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[7] + a[3]*b[6] + a[4]*b[13] + a[5]*b[12];
//		dest[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[8] - a[3]*b[9] + a[4]*b[14] - a[5]*b[15];
//		dest[3] = a[0]*b[3] + a[1]*b[2] + a[2]*b[9] + a[3]*b[8] + a[4]*b[15] + a[5]*b[14];
//		dest[4] = a[0]*b[4] - a[1]*b[5] + a[2]*b[10] - a[3]*b[11] + a[4]*b[16] - a[5]*b[17];
//		dest[5] = a[0]*b[5] + a[1]*b[4] + a[2]*b[11] + a[3]*b[10] + a[4]*b[17] + a[5]*b[16];
//		dest[6] = a[6]*b[0] - a[7]*b[1] + a[8]*b[6] - a[9]*b[7] + a[10]*b[12] - a[11]*b[13];
//		dest[7] = a[6]*b[1] + a[7]*b[0] + a[8]*b[7] + a[9]*b[6] + a[10]*b[13] + a[11]*b[12];
//		dest[8] = a[6]*b[2] - a[7]*b[3] + a[8]*b[8] - a[9]*b[9] + a[10]*b[14] - a[11]*b[15];
//		dest[9] = a[6]*b[3] + a[7]*b[2] + a[8]*b[9] + a[9]*b[8] + a[10]*b[15] + a[11]*b[14];
//		dest[10] = a[6]*b[4] - a[7]*b[5] + a[8]*b[10] - a[9]*b[11] + a[10]*b[16] - a[11]*b[17];
//		dest[11] = a[6]*b[5] + a[7]*b[4] + a[8]*b[11] + a[9]*b[10] + a[10]*b[17] + a[11]*b[16];
//		dest[12] = a[12]*b[0] - a[13]*b[1] + a[14]*b[6] - a[15]*b[7] + a[16]*b[12] - a[17]*b[13];
//		dest[13] = a[12]*b[1] + a[13]*b[0] + a[14]*b[7] + a[15]*b[6] + a[16]*b[13] + a[17]*b[12];
//		dest[14] = a[12]*b[2] - a[13]*b[3] + a[14]*b[8] - a[15]*b[9] + a[16]*b[14] - a[17]*b[15];
//		dest[15] = a[12]*b[3] + a[13]*b[2] + a[14]*b[9] + a[15]*b[8] + a[16]*b[15] + a[17]*b[14];
//		dest[16] = a[12]*b[4] - a[13]*b[5] + a[14]*b[10] - a[15]*b[11] + a[16]*b[16] - a[17]*b[17];
//		dest[17] = a[12]*b[5] + a[13]*b[4] + a[14]*b[11] + a[15]*b[10] + a[16]*b[17] + a[17]*b[16];
//	}
//
//	static CUDA_HOST_DEVICE void inline flipSign( TYPE& a )
//	{
//		a = -a;
//	}
//
//	static CUDA_HOST_DEVICE void inline swap( TYPE& a, TYPE& b )
//	{
//		TYPE temp = b;
//		b = a;
//		a = temp;
//	}
//
//	static CUDA_HOST_DEVICE void inline swapAndFlipSign( TYPE& a, TYPE& b )
//	{
//		TYPE temp = b;
//		b = -a;
//		a = -temp;
//	}
//
//	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
//	{
//		flipSign( store[1] );
//		swap( store[2], store[6] );
//		swapAndFlipSign( store[3], store[7] );
//		swap( store[4], store[12] );
//		swapAndFlipSign( store[5], store[13] );
//		flipSign( store[9] );
//		swap( store[10], store[14] );
//		swapAndFlipSign( store[11], store[15] );
//		flipSign( store[17] );
//	}

	static CUDA_HOST_DEVICE inline lat_group_index_t getIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return i*NC+j;
	}

	static CUDA_HOST_DEVICE typename Real4<REALTYPE>::VECTORTYPE inline getSU2Subgroup( TYPE store[SIZE], lat_group_index_t i, lat_group_index_t j )
	{
		typename Real4<REALTYPE>::VECTORTYPE result;

		TYPE temp;

		temp = store[getIndex(i,i)];
		result.x = temp.x;
		result.w = temp.y;

		temp = store[getIndex(j,i)];
		result.x += temp.x;
		result.w -= temp.y;

		temp = store[getIndex(i,j)];
		result.z = temp.x;
		result.y = temp.y;

		temp = store[getIndex(j,j)];
		result.z -= temp.x;
		result.y += temp.y;

		return result;
	}


	static CUDA_HOST_DEVICE void inline rightSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		for( lat_group_index_t k = 0; k < 3; k++ )
		{
			TYPE KI = TYPE( mat.x, mat.y ) * store[getIndex(k,i)];
			KI += TYPE( -mat.z, mat.w ) * store[getIndex(k,j)];

			TYPE KJ = TYPE( mat.z, mat.w ) * store[getIndex(k,i)];
			KJ += TYPE( mat.x, -mat.y ) * store[getIndex(k,j)];

			store[getIndex(k,i)] = KI;
			store[getIndex(k,j)] = KJ;
		}
	}

	static CUDA_HOST_DEVICE void inline leftSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		for( lat_group_index_t k = 0; k < 3; k++ )
		{
			TYPE IK = TYPE( mat.x, mat.y ) * store[getIndex(i,k)];
			IK += TYPE( mat.z, mat.w ) * store[getIndex(j,k)];

			TYPE JK = TYPE( -mat.z, mat.w ) * store[getIndex(i,k)];
			JK += TYPE( mat.x, -mat.y ) * store[getIndex(j,k)];

			store[getIndex(i,k)] = IK;
			store[getIndex(j,k)] = JK;
		}
	}
};


} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
