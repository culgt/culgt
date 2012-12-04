/*
 * RandomUpdate.hxx
 *
 *  Created on: May 10, 2012
 *      Author: vogt
 */

#ifndef RANDOMUPDATE_HXX_
#define RANDOMUPDATE_HXX_

#include "../../../util/datatype/datatypes.h"
#include "../../../util/rng/PhiloxWrapper.hxx"

class RandomUpdate
{
public:
	__device__ inline RandomUpdate();
	__device__ inline RandomUpdate( PhiloxWrapper *rng );
	__device__ inline void calculateUpdate( volatile Real (&shA)[128], short id );
private:
	PhiloxWrapper *rng;
};

__device__ RandomUpdate::RandomUpdate()
{
}

__device__ RandomUpdate::RandomUpdate( PhiloxWrapper *rng ) : rng(rng)
{
}

__device__ void RandomUpdate::calculateUpdate( volatile Real (&shA)[128], short id )
{
	Real alpha, phi, cos_theta, sin_theta, sin_alpha;
	alpha = rng->rand();
	phi = 2.0 * rng->rand();
	cos_theta = 2.0 * rng->rand() - 1.0;
	sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	sin_alpha = sinpi(alpha);
	shA[id] = cospi(alpha);
	shA[id+NSB] = sin_alpha * sin_theta * cospi(phi);
	shA[id+NSB2] = sin_alpha * sin_theta * sinpi(phi);
	shA[id+NSB3] = sin_alpha * cos_theta;
}
#endif /* ORUPDATE_HXX_ */
