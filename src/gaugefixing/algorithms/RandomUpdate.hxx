/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 */

#ifndef RANDOMUPDATE_HXX_
#define RANDOMUPDATE_HXX_

#include "../../lattice/datatype/datatypes.h"
#include "../../lattice/rng/PhiloxWrapper.hxx"

class RandomUpdate
{
public:
	__device__ inline RandomUpdate();
	__device__ inline RandomUpdate( PhiloxWrapper *rng );
	__device__ inline void calculateUpdate( volatile Real (&shA)[4*NSB], short id );
private:
	PhiloxWrapper *rng;
};

__device__ RandomUpdate::RandomUpdate()
{
}

__device__ RandomUpdate::RandomUpdate( PhiloxWrapper *rng ) : rng(rng)
{
}

__device__ void RandomUpdate::calculateUpdate( volatile Real (&shA)[4*NSB], short id )
{
	Real alpha, phi, cos_theta, sin_theta, sin_alpha;
	alpha = rng->rand();
	phi = 2.0 * rng->rand();
	cos_theta = 2.0 * rng->rand() - 1.0;
	sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	sin_alpha = sinpi(alpha);
	shA[id] = cospi(alpha);
	shA[id+NSB] = sin_alpha * sin_theta * cospi(phi);
	shA[id+2*NSB] = sin_alpha * sin_theta * sinpi(phi);
	shA[id+3*NSB] = sin_alpha * cos_theta;
}
#endif /* ORUPDATE_HXX_ */
