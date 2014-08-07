/**
 * Flops: 14
 *
 */

#ifndef MICROUPDATE_H_
#define MICROUPDATE_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

template<typename T, typename InternalT=T>  class MicroUpdate
{
public:
	__device__ inline MicroUpdate();
	__device__ inline void calculateUpdate( T* shA, const short& id, const int NSB );
	const static int Flops = 14;
private:
};

template<typename T, typename InternalT> __device__ MicroUpdate<T,InternalT>::MicroUpdate()
{
}

template<typename T, typename InternalT> __device__ void MicroUpdate<T,InternalT>::calculateUpdate( T* shA, const short& id, const int NSB )
{
	InternalT ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	InternalT a0_sq = shA[id]*shA[id];

	InternalT b=2.*shA[id]/(a0_sq+ai_sq);

	shA[id]=(a0_sq-ai_sq)/(a0_sq+ai_sq);
	shA[id+NSB]*=b;
	shA[id+2*NSB]*=b;
	shA[id+3*NSB]*=b;
}

}

#endif
