/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU2REALFULLIGNORESECONDLINE_H_
#define SU2REALFULLIGNORESECONDLINE_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"
#include "../../cuLGT1legacy/Complex.hxx"

namespace culgt
{

template<typename T> class SU2RealFullIgnoreSecondLine
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 8;
	typedef T TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE/2; i++ )
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
	}

	static CUDA_HOST_DEVICE Complex<REALTYPE> inline trace( TYPE store[SIZE] )
	{
		return Complex<REALTYPE>( 2.*store[0], 2.*store[1] );
	}

	static CUDA_HOST_DEVICE inline REALTYPE reDet( TYPE store[SIZE] )
	{
		return store[0]*store[0] + store[1]*store[1] + store[2]*store[2] + store[3]*store[3];
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return 2.*store[0];
	}

	static CUDA_HOST_DEVICE void inline multAssignScalar( TYPE dest[SIZE], const REALTYPE scalar )
	{
		for( int i = 0; i < SIZE/2; i++ )
		{
			dest[i] *= scalar;
		}
	}

	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		TYPE a[SIZE];
		for( int i = 0; i < SIZE/2; i++ )
		{
			a[i] = dest[i];
		}

		dest[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
		dest[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
		dest[2] = a[0]*b[2] + a[2]*b[0] - a[1]*b[3] + a[3]*b[1];
		dest[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
	}

	static CUDA_HOST_DEVICE void inline addAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] += b[0];
		dest[1] += b[1];
		dest[2] += b[2];
		dest[3] += b[3];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0] -= b[0];
		dest[1] -= b[1];
		dest[2] -= b[2];
		dest[3] -= b[3];
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const Complex<REALTYPE> b )
	{
		dest[0] -= b.x;
		dest[1] -= b.y;
	}

	static CUDA_HOST_DEVICE REALTYPE inline normFrobeniusSquared( TYPE store[SIZE] )
	{
		return store[0]*store[0]+store[1]*store[1]+store[2]*store[2]+store[3]*store[3];
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
		flipSign( store[2] );
		flipSign( store[3] );
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


} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
