/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU3VECTOR2_H_
#define SU3VECTOR2_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"

namespace culgt
{



/*
 *
 */
template<typename T> class SU3Vector2
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = 6;
	typedef typename Real2<T>::VECTORTYPE TYPE;
	typedef T REALTYPE;

	static CUDA_HOST_DEVICE void inline zero( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[0].x = 0;
			store[0].y = 0;
		}
	}

	static CUDA_HOST_DEVICE void inline identity( TYPE store[SIZE] )
	{
		zero( store );
		store[0].x = 1;
		store[4].x = 1;
	}
};

} /* namespace culgt */
#endif /* SU[0].wVECTOR4_H_ */
