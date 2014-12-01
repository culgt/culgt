/**
 *
 */

#ifndef RANDOMUPDATE_H_
#define RANDOMUPDATE_H_

#include "../../common/culgt_typedefs.h"
#include "../../util/rng/PhiloxWrapper.h"

namespace culgt
{

template<typename T, typename InternalT=T> class RandomUpdate
{
public:
	__device__ inline RandomUpdate( PhiloxWrapper<InternalT>& rng ) : rng(rng){};
	__device__ inline void calculateUpdate( T* shA, const short& id, const int NSB );
	const static int Flops = -1;
private:
	PhiloxWrapper<InternalT> rng;
};


template<typename T, typename InternalT> __device__ void RandomUpdate<T,InternalT>::calculateUpdate( T* shA, const short& id, const int NSB )
{
	InternalT alpha, phi, cos_theta, sin_theta, sin_alpha;
	alpha = rng.rand();
	phi = 2.0 * rng.rand();
	cos_theta = 2.0 * rng.rand() - 1.0;
	sin_theta = ::sqrt(1.0 - cos_theta * cos_theta);
	sin_alpha = sinpi(alpha);
	shA[id] = cospi(alpha);
	shA[id+NSB] = sin_alpha * sin_theta * cospi(phi);
	shA[id+2*NSB] = sin_alpha * sin_theta * sinpi(phi);
	shA[id+3*NSB] = sin_alpha * cos_theta;
}

}
#endif /* ORUPDATE_HXX_ */
