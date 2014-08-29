/**
 *
 * Flops: 22
 *
 */

#ifndef ORUPDATE_H_
#define ORUPDATE_H_

#include "../../common/culgt_typedefs.h"

template<typename T, typename InternalT=T> class OrUpdate
{
public:
	__device__ inline OrUpdate( T orParameter ) : orParameter(orParameter){};
	__device__ inline void calculateUpdate( T* shA, const short& id, const int NSB );
	__device__ inline void calculateUpdate( T* shA, typename Real4<T>::VECTORTYPE& q, const short& id, const int NSB );
	const static int Flops = 22;
private:
	T orParameter;
};


template<typename T, typename InternalT> __device__ void OrUpdate<T,InternalT>::calculateUpdate( T* shA, const short& id, const int NSB )
{
	InternalT ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	InternalT a0_sq = shA[id]*shA[id];

	InternalT b=((InternalT)orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
	InternalT c=rsqrt(a0_sq+b*b*ai_sq);

	shA[id]*=c;
	shA[id+NSB]*=b*c;
	shA[id+2*NSB]*=b*c;
	shA[id+3*NSB]*=b*c;

	// 22 flops
}

template<typename T, typename InternalT> __device__ void OrUpdate<T,InternalT>::calculateUpdate( T* shA, typename Real4<T>::VECTORTYPE& q, const short& id, const int NSB )
{
	q.x = shA[id];
	q.y = shA[id+NSB];
	q.z = shA[id+2*NSB];
	q.w = shA[id+3*NSB];

	InternalT ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	InternalT a0_sq = shA[id]*shA[id];

	InternalT b=((InternalT)orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
	InternalT c=rsqrt(a0_sq+b*b*ai_sq);

	q.x*=c;
	q.y*=b*c;
	q.z*=b*c;
	q.w*=b*c;

	// 22 flops
}

#endif /* ORUPDATE_HXX_ */
