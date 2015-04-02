/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SUNREALFULL_H_
#define SUNREALFULL_H_

#include "common/culgt_typedefs.h"
#include "cudacommon/cuda_host_device.h"
#include "math/Complex.h"

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

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
	}
};

template<typename T> class SUNRealFull<2,T>
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 8;
	typedef T TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
	}

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		store[0] = 1.;
		store[1] = 0.;
		store[2] = 0.;
		store[3] = 0.;
		store[4] = 0.;
		store[5] = 0.;
		store[6] = 1.;
		store[7] = 0.;
	}

	static CUDA_HOST_DEVICE Complex<REALTYPE> inline trace( TYPE store[SIZE] )
	{
		return Complex<REALTYPE>( store[0]+store[6], store[1]+store[7] );
	}

	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
	{
		return store[0]*store[6] - store[2]*store[4] - store[1]*store[7] + store[3]*store[5];
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return store[0]+store[6];
	}

	static CUDA_HOST_DEVICE void inline multAssignScalar( TYPE dest[SIZE], const REALTYPE scalar )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			dest[i] *= scalar;
		}
	}

	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		TYPE a[SIZE];
		for( int i = 0; i < SIZE; i++ )
		{
			a[i] = dest[i];
		}

		dest[0] = a[0]*b[0] - a[1]*b[1] + a[2]*b[4] - a[3]*b[5];
		dest[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[5] + a[3]*b[4];
		dest[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[6] - a[3]*b[7];
		dest[3] = a[0]*b[3] + a[1]*b[2] + a[2]*b[7] + a[3]*b[6];
		dest[4] = a[4]*b[0] - a[5]*b[1] + a[6]*b[4] - a[7]*b[5];
		dest[5] = a[4]*b[1] + a[5]*b[0] + a[6]*b[5] + a[7]*b[4];
		dest[6] = a[4]*b[2] - a[5]*b[3] + a[6]*b[6] - a[7]*b[7];
		dest[7] = a[4]*b[3] + a[5]*b[2] + a[6]*b[7] + a[7]*b[6];
	}

	static CUDA_HOST_DEVICE void inline addAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] += b[0];
		dest[1] += b[1];
		dest[2] += b[2];
		dest[3] += b[3];
		dest[4] += b[4];
		dest[5] += b[5];
		dest[6] += b[6];
		dest[7] += b[7];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] -= b[0];
		dest[1] -= b[1];
		dest[2] -= b[2];
		dest[3] -= b[3];
		dest[4] -= b[4];
		dest[5] -= b[5];
		dest[6] -= b[6];
		dest[7] -= b[7];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const Complex<REALTYPE> b )
	{
		dest[0] -= b.x;
		dest[1] -= b.y;
		dest[6] -= b.x;
		dest[7] -= b.y;
	}

	static CUDA_HOST_DEVICE REALTYPE inline normFrobeniusSquared( TYPE store[SIZE] )
	{
		return store[0]*store[0]+store[1]*store[1]+store[2]*store[2]+store[3]*store[3] +store[4]*store[4] +store[5]*store[5]+store[6]*store[6]+store[7]*store[7];
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
		swap( store[2], store[4] );
		swapAndFlipSign( store[3], store[5] );
		flipSign( store[7] );
	}

	static CUDA_HOST_DEVICE void inline reproject( TYPE store[SIZE] )
	{
		T fac = 1./::sqrt( reDet( store ) );

		store[0] *= fac;
		store[1] *= fac;
		store[2] *= fac;
		store[3] *= fac;

		store[4] = -store[2];
		store[5] *= store[3];
		store[6] *= store[0];
		store[7] *= -store[1];
	}

	/**
	 * Ignores second line
	 */
	static const int FlopsGetSU2Subgroup = 0;
	static CUDA_HOST_DEVICE typename Real4<REALTYPE>::VECTORTYPE inline getSU2Subgroup( TYPE store[SIZE], lat_group_index_t i, lat_group_index_t j )
	{
		typename Real4<REALTYPE>::VECTORTYPE result;
		result.x = store[0];
		result.y = store[1];
		result.z = store[2];
		result.w = store[3];
		return result;
	}

	static const int FlopsSubgroupMult = 4*7;
	static CUDA_HOST_DEVICE void inline rightSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		store[6] = store[0]*mat.x - store[3]*mat.w - store[1]*mat.y - store[2]*mat.z;
		store[7] = -(store[2]*mat.w + store[1]*mat.x + store[0]*mat.y - store[3]*mat.z);
		store[4] = -(store[2]*mat.x - store[1]*mat.w + store[3]*mat.y + store[0]*mat.z);
		store[5] = store[0]*mat.w + store[3]*mat.x - store[2]*mat.y + store[1]*mat.z;
		store[2] = -store[4];
		store[3] = store[5];
		store[0] = store[6];
		store[1] = -store[7];
	}

	static CUDA_HOST_DEVICE void inline leftSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		store[6] = store[0]*mat.x - store[3]*mat.w - store[1]*mat.y - store[2]*mat.z;
		store[7] = -(store[1]*mat.x - store[2]*mat.w + store[0]*mat.y + store[3]*mat.z);
		store[4] = -(store[1]*mat.w + store[2]*mat.x - store[3]*mat.y + store[0]*mat.z);
		store[5] = store[0]*mat.w + store[3]*mat.x + store[2]*mat.y - store[1]*mat.z;
		store[2] = -store[4];
		store[3] = store[5];
		store[0] = store[6];
		store[1] = -store[7];
	}

};

template<typename T> class SUNRealFull<3,T>
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = 2*(3*3);
	typedef T TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
	}

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

	static CUDA_HOST_DEVICE REALTYPE inline normFrobeniusSquared( TYPE store[SIZE] )
	{
		REALTYPE result = 0;
		for( int i = 0; i < SIZE; i++ )
		{
			result += store[i]*store[i];
		}
		return result;
	}

	static CUDA_HOST_DEVICE Complex<REALTYPE> inline trace( TYPE store[SIZE] )
	{
		return Complex<REALTYPE>(store[0]+store[8]+store[16], store[1]+store[9]+store[17]);
	}

	static CUDA_HOST_DEVICE void inline multAssignScalar( TYPE dest[SIZE], const REALTYPE scalar )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			dest[i] *= scalar;
		}
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

	static CUDA_HOST_DEVICE void inline addAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] += b[0];
		dest[1] += b[1];
		dest[2] += b[2];
		dest[3] += b[3];
		dest[4] += b[4];
		dest[5] += b[5];
		dest[6] += b[6];
		dest[7] += b[7];
		dest[8] += b[8];
		dest[9] += b[9];
		dest[10] += b[10];
		dest[11] += b[11];
		dest[12] += b[12];
		dest[13] += b[13];
		dest[14] += b[14];
		dest[15] += b[15];
		dest[16] += b[16];
		dest[17] += b[17];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] -= b[0];
		dest[1] -= b[1];
		dest[2] -= b[2];
		dest[3] -= b[3];
		dest[4] -= b[4];
		dest[5] -= b[5];
		dest[6] -= b[6];
		dest[7] -= b[7];
		dest[8] -= b[8];
		dest[9] -= b[9];
		dest[10] -= b[10];
		dest[11] -= b[11];
		dest[12] -= b[12];
		dest[13] -= b[13];
		dest[14] -= b[14];
		dest[15] -= b[15];
		dest[16] -= b[16];
		dest[17] -= b[17];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const Complex<REALTYPE> b )
	{
		dest[0] -= b.x;
		dest[1] -= b.y;
		dest[8] -= b.x;
		dest[9] -= b.y;
		dest[16] -= b.x;
		dest[17] -= b.y;
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

	static CUDA_HOST_DEVICE lat_group_index_t inline getRealIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return (3*i+j)*2;
	}
	static CUDA_HOST_DEVICE lat_group_index_t inline getImagIndex( lat_group_index_t i, lat_group_index_t j )
	{
		return (3*i+j)*2+1;
	}

	static const int FlopsGetSU2Subgroup = 4;
	static CUDA_HOST_DEVICE typename Real4<REALTYPE>::VECTORTYPE inline getSU2Subgroup( TYPE store[SIZE], lat_group_index_t i, lat_group_index_t j )
	{
		typename Real4<REALTYPE>::VECTORTYPE result;
		result.x = store[getRealIndex(i,i)]+store[getRealIndex(j,j)];
		result.w = store[getImagIndex(i,i)]-store[getImagIndex(j,j)];
		result.z = store[getRealIndex(i,j)]-store[getRealIndex(j,i)];
		result.y = store[getImagIndex(i,j)]+store[getImagIndex(j,i)];
		return result;
	}

	static const int FlopsSubgroupMult = 3*4*7;
	static CUDA_HOST_DEVICE void inline rightSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		for( lat_group_index_t k = 0; k < 3; k++ )
		{
			TYPE KIr = mat.x*store[getRealIndex(k,i)] - mat.w*store[getImagIndex(k,i)] - mat.z*store[getRealIndex(k,j)] - mat.y*store[getImagIndex(k,j)];
			TYPE KIi = mat.x*store[getImagIndex(k,i)] + mat.w*store[getRealIndex(k,i)] - mat.z*store[getImagIndex(k,j)] + mat.y*store[getRealIndex(k,j)];

			TYPE KJr = +mat.z*store[getRealIndex(k,i)] - mat.y*store[getImagIndex(k,i)] + mat.x*store[getRealIndex(k,j)] + mat.w*store[getImagIndex(k,j)];
			TYPE KJi = +mat.z*store[getImagIndex(k,i)] + mat.y*store[getRealIndex(k,i)] + mat.x*store[getImagIndex(k,j)] - mat.w*store[getRealIndex(k,j)];

			store[getRealIndex(k,i)] = KIr;
			store[getImagIndex(k,i)] = KIi;

			store[getRealIndex(k,j)] = KJr;
			store[getImagIndex(k,j)] = KJi;
		}
	}

	static CUDA_HOST_DEVICE void inline leftSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		for( lat_group_index_t k = 0; k < 3; k++ )
		{
			TYPE IKr = mat.x*store[getRealIndex(i,k)] - mat.w*store[getImagIndex(i,k)] + mat.z*store[getRealIndex(j,k)] - mat.y*store[getImagIndex(j,k)];
			TYPE IKi = mat.x*store[getImagIndex(i,k)] + mat.w*store[getRealIndex(i,k)] + mat.z*store[getImagIndex(j,k)] + mat.y*store[getRealIndex(j,k)];

			TYPE JKr = -mat.z*store[getRealIndex(i,k)] - mat.y*store[getImagIndex(i,k)] + mat.x*store[getRealIndex(j,k)] + mat.w*store[getImagIndex(j,k)];
			TYPE JKi = -mat.z*store[getImagIndex(i,k)] + mat.y*store[getRealIndex(i,k)] + mat.x*store[getImagIndex(j,k)] - mat.w*store[getRealIndex(j,k)];

			store[getRealIndex(i,k)] = IKr;
			store[getImagIndex(i,k)] = IKi;

			store[getRealIndex(j,k)] = JKr;
			store[getImagIndex(j,k)] = JKi;
		}
	}

	static CUDA_HOST_DEVICE void inline normalizeRow( int row, TYPE store[SIZE] )
	{
		int baseIndex = row*6;
		T norm = store[baseIndex]*store[baseIndex]+store[baseIndex+1]*store[baseIndex+1]+store[baseIndex+2]*store[baseIndex+2]+store[baseIndex+3]*store[baseIndex+3]+store[baseIndex+4]*store[baseIndex+4]+store[baseIndex+5]*store[baseIndex+5];
		norm = 1./::sqrt( norm ); // TODO replace by rsqrt and make compatible with non-cuda compiler

		store[baseIndex+0] *= norm;
		store[baseIndex+1] *= norm;
		store[baseIndex+2] *= norm;
		store[baseIndex+3] *= norm;
		store[baseIndex+4] *= norm;
		store[baseIndex+5] *= norm;
	}

	static CUDA_HOST_DEVICE void inline reproject( TYPE store[SIZE] )
	{
		// normalize first row
		normalizeRow( 0, store );
//		T norm = store[0]*store[0]+store[1]*store[1]+store[2]*store[2]+store[3]*store[3]+store[4]*store[4]+store[5]*store[5];
//		norm = 1./::sqrt( norm ); // TODO replace by rsqrt and make compatible with non-cuda compiler
//
//		store[0] *= norm;
//		store[1] *= norm;
//		store[2] *= norm;
//		store[3] *= norm;
//		store[4] *= norm;
//		store[5] *= norm;

		T a6 = store[6] - store[0]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) + store[1]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);
		T a7 = store[7] - store[1]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) - store[0]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);
		T a8 = store[8] - store[2]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) + store[3]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);
		T a9 = store[9] - store[3]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) - store[2]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);
		T a10 = store[10] - store[4]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) + store[5]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);
		T a11 = store[11] - store[5]*(store[0]*store[6] + store[1]*store[7] + store[2]*store[8] + store[3]*store[9] + store[4]*store[10] + store[5]*store[11]) - store[4]*(store[0]*store[7] - store[1]*store[6] + store[2]*store[9] - store[3]*store[8] + store[4]*store[11] - store[5]*store[10]);

		store[6] = a6;
		store[7] = a7;
		store[8] = a8;
		store[9] = a9;
		store[10] = a10;
		store[11] = a11;

		normalizeRow( 1, store );

		store[12] = store[2]*store[10]-store[3]*store[11] - store[8]*store[4]+store[9]*store[5];
		store[13] = -store[2]*store[11]-store[3]*store[10] + store[8]*store[5]+store[9]*store[4];
		store[14] = store[4]*store[6]-store[5]*store[7] - store[10]*store[0]+store[11]*store[1] ;
		store[15] = -store[4]*store[7]-store[5]*store[6] + store[10]*store[1]+store[11]*store[0];
		store[16] = store[0]*store[8]-store[1]*store[9] - store[2]*store[6]+store[3]*store[7];
		store[17] = -store[0]*store[9]-store[1]*store[8] + store[2]*store[7]+store[3]*store[6];

		normalizeRow( 2, store );

//		T norm = a12*a12+a13*a13+a14*a14+a15*a15+a16*a16+a17*a17;
//		norm = 1./::sqrt( norm ); // TODO replace by rsqrt and make compatible with non-cuda compiler
//
//		l1.set(12, a12*norm );
//		l1.set(13, a13*norm );
//		l1.set(14, a14*norm );
//		l1.set(15, a15*norm );
//		l1.set(16, a16*norm );
//		l1.set(17, a17*norm );

	}
};


} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
