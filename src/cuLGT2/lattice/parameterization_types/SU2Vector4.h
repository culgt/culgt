/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU2VECTOR4_H_
#define SU2VECTOR4_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename T> struct Vector4
{
};

template<> struct Vector4<float>
{
	typedef float4 VECTORTYPE;
};

template<> struct Vector4<double>
{
	typedef double4 VECTORTYPE;
};

/*
 *
 */
template<typename T> class SU2Vector4
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 1;
	typedef typename Vector4<T>::VECTORTYPE TYPE;
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

	static CUDA_HOST_DEVICE REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return 2.*store[0].x;
	}

	static CUDA_HOST_DEVICE REALTYPE inline reDet( TYPE store[SIZE] )
	{
		return store[0].x*store[0].x+store[0].y*store[0].y+store[0].z*store[0].z+store[0].w*store[0].w;
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

	static CUDA_HOST_DEVICE void inline hermitian( TYPE store[SIZE] )
	{
		store[0].y = -store[0].y;
		store[0].z = -store[0].z;
		store[0].w = -store[0].w;
	}
};

} /* namespace culgt */
#endif /* SU[0].wVECTOR4_H_ */
