/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU2VECTOR4_H_
#define SU2VECTOR4_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"
#include "../../cuLGT1legacy/Complex.hxx"

namespace culgt
{



/*
 *
 */
template<typename T> class SU2Vector4
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 1;
	typedef typename Real4<T>::VECTORTYPE TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		store[0].x = 0;
		store[0].y = 0;
		store[0].z = 0;
		store[0].w = 0;
	}

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		store[0].x = 1;
		store[0].y = 0;
		store[0].z = 0;
		store[0].w = 0;
	}

	static CUDA_HOST_DEVICE Complex<REALTYPE> inline trace( TYPE store[SIZE] )
	{
		return Complex<REALTYPE>( 2.*store[0].x, 0. );
	}

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return 2.*store[0].x;
	}

	static CUDA_HOST_DEVICE REALTYPE inline reDet( TYPE store[SIZE] )
	{
		return store[0].x*store[0].x+store[0].y*store[0].y+store[0].z*store[0].z+store[0].w*store[0].w;
	}

	static CUDA_HOST_DEVICE void inline multAssignScalar( TYPE dest[SIZE], const REALTYPE scalar )
	{
		dest[0].x *= scalar;
		dest[0].y *= scalar;
		dest[0].z *= scalar;
		dest[0].w *= scalar;
	}

	static CUDA_HOST_DEVICE void inline multAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		TYPE a[SIZE];

		a[0] = dest[0];

		dest[0].x = a[0].x*b[0].x - a[0].y*b[0].y - a[0].z*b[0].z - a[0].w*b[0].w;
		dest[0].y = a[0].x*b[0].y + a[0].y*b[0].x + a[0].z*b[0].w - a[0].w*b[0].z;
		dest[0].z = a[0].x*b[0].z + a[0].z*b[0].x - a[0].y*b[0].w + a[0].w*b[0].y;
		dest[0].w = a[0].x*b[0].w + a[0].y*b[0].z - a[0].z*b[0].y + a[0].w*b[0].x;
	}

	static CUDA_HOST_DEVICE void inline addAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0].x += b[0].x;
		dest[0].y += b[0].y;
		dest[0].z += b[0].z;
		dest[0].w += b[0].w;
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const TYPE b[SIZE] )
	{
		dest[0].x -= b[0].x;
		dest[0].y -= b[0].y;
		dest[0].z -= b[0].z;
		dest[0].w -= b[0].w;
	}

	static CUDA_HOST_DEVICE void inline subtractAssign( TYPE dest[SIZE], const Complex<REALTYPE> b )
	{
		dest[0].x -= b.x;
		dest[0].y -= b.y;
	}

	static CUDA_HOST_DEVICE REALTYPE inline normFrobeniusSquared( TYPE store[SIZE] )
	{
		return 2.*(store[0].x*store[0].x+store[0].y*store[0].y+store[0].z*store[0].z+store[0].w*store[0].w);
	}

	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
	{
		store[0].y = -store[0].y;
		store[0].z = -store[0].z;
		store[0].w = -store[0].w;
	}



	/**
	 * The following subgroup methods are defined for convenience.
	 */
	static const int FlopsGetSU2Subgroup = 0;
	static CUDA_HOST_DEVICE typename Real4<REALTYPE>::VECTORTYPE inline getSU2Subgroup( TYPE store[SIZE], lat_group_index_t i, lat_group_index_t j )
	{
		return store[0];
	}

	static const int FlopsSubgroupMult = 4*7;
	static CUDA_HOST_DEVICE void inline rightSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		TYPE a[SIZE];
		a[0] = store[0];

		store[0].x = a[0].x*mat.x - a[0].y*mat.y - a[0].z*mat.z - a[0].w*mat.w;
		store[0].y = a[0].x*mat.y + a[0].y*mat.x + a[0].z*mat.w - a[0].w*mat.z;
		store[0].z = a[0].x*mat.z + a[0].z*mat.x - a[0].y*mat.w + a[0].w*mat.y;
		store[0].w = a[0].x*mat.w + a[0].y*mat.z - a[0].z*mat.y + a[0].w*mat.x;
	}

	static CUDA_HOST_DEVICE void inline leftSubgroupMult( TYPE store[SIZE], const typename Real4<REALTYPE>::VECTORTYPE& mat, lat_group_index_t i, lat_group_index_t j )
	{
		TYPE a[SIZE];
		a[0] = store[0];

		store[0].x = a[0].x*mat.x - a[0].y*mat.y - a[0].z*mat.z - a[0].w*mat.w;
		store[0].y = a[0].y*mat.x + a[0].x*mat.y + a[0].w*mat.z - a[0].z*mat.w;
		store[0].z = a[0].x*mat.z + a[0].z*mat.x - a[0].w*mat.y + a[0].y*mat.w;
		store[0].w = a[0].x*mat.w + a[0].w*mat.x - a[0].y*mat.z + a[0].z*mat.y;
	}

};

} /* namespace culgt */
#endif /* SU[0].wVECTOR4_H_ */
