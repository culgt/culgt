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

template<typename T> class SUNComplexFull<2,T>
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 4;
	typedef T TYPE;
	typedef T REALTYPE;

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
};

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

	static CUDA_HOST_DEVICE REALTYPE inline reDet( TYPE store[SIZE] )
	{
		Complex<T> result(0.,0.);

		result += store[0] * store[4] * store[8];
		result += store[1] * store[5] * store[6];
		result += store[2] * store[3] * store[7];
		result -= store[2] * store[4] * store[6];
		result -= store[0] * store[5] * store[7];
		result -= store[1] * store[3] * store[8];

		return result.x;
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
	}

	static CUDA_HOST_DEVICE void inline multAssignScalarComplex( TYPE dest[SIZE], const Complex<REALTYPE> scalar )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			dest[i] *= scalar;
		}
	}

//	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
//	{
//		return store[0]*store[8]*store[16] - store[0]*store[10]*store[14] - store[2]*store[6]*store[16] + store[2]*store[10]*store[12] + store[4]*store[6]*store[14] - store[4]*store[8]*store[12] - store[0]*store[9]*store[17] + store[0]*store[11]*store[15] - store[1]*store[8]*store[17] - store[1]*store[9]*store[16] + store[1]*store[10]*store[15] + store[1]*store[11]*store[14] + store[2]*store[7]*store[17] - store[2]*store[11]*store[13] + store[3]*store[6]*store[17] + store[3]*store[7]*store[16] - store[3]*store[10]*store[13] - store[3]*store[11]*store[12] - store[4]*store[7]*store[15] + store[4]*store[9]*store[13] - store[5]*store[6]*store[15] - store[5]*store[7]*store[14] + store[5]*store[8]*store[13] + store[5]*store[9]*store[12];
//	}
//
	static CUDA_HOST_DEVICE Complex<REALTYPE> inline trace( TYPE store[SIZE] )
	{
		return store[0]+store[4]+store[8];
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return store[0].x+store[4].x+store[8].x;
	}

	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		TYPE a[SIZE];
		for( int i = 0; i < SIZE; i++ )
		{
			a[i] = dest[i];
		}

		for( int i = 0; i < 3; i++ )
		{
			for( int j = 0; j < 3; j++ )
			{
				dest[i*3+j] = TYPE(0,0);
				for( int k = 0; k < 3; k++ )
				{
					dest[i*3+j] += a[i*3+k]*b[k*3+j];
				}
			}
		}
	}


	static CUDA_HOST_DEVICE void inline swapAndConjugate( TYPE& a, TYPE& b )
	{
		TYPE temp = b;
		b = a.conj();
		a = temp.conj();
	}

	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
	{
		store[0] = store[0].conj();
		swapAndConjugate( store[1], store[3] );
		swapAndConjugate( store[2], store[6] );
		store[4] = store[4].conj();
		swapAndConjugate( store[5], store[7] );
		store[8] = store[8].conj();
	}

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
