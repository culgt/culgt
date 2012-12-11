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
 *
 * Flops: 22
 *
 */

#ifndef ORUPDATE_HXX_
#define ORUPDATE_HXX_

#include "../GlobalConstants.h"
#include "../../lattice/datatype/datatypes.h"

class OrUpdate
{
public:
	__device__ inline OrUpdate();
	__device__ inline OrUpdate( float param );
	__device__ inline void calculateUpdate( volatile Real (&shA)[4*NSB], short id );
	__device__ inline void setParameter( float param );
	__device__ inline float getParameter();
private:
	float orParameter;
};

__device__ OrUpdate::OrUpdate()
{
}

__device__ OrUpdate::OrUpdate( float param ) : orParameter(param)
{
}

__device__ void OrUpdate::calculateUpdate( volatile Real (&shA)[4*NSB], short id )
{
#ifdef USE_DP_ORUPDATE
	double ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	double a0_sq = shA[id]*shA[id];

	double b=((double)orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
	double c=rsqrt(a0_sq+b*b*ai_sq); // even in DP calculation the rsqrt has no effect on the stability of GFF

	shA[id]*=c;
	shA[id+NSB]*=b*c;
	shA[id+2*NSB]*=b*c;
	shA[id+3*NSB]*=b*c;
#else
	Real ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	Real a0_sq = shA[id]*shA[id];

	Real b=(orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
	Real c=rsqrt(a0_sq+b*b*ai_sq);

	shA[id]*=c;
	shA[id+NSB]*=b*c;
	shA[id+2*NSB]*=b*c;
	shA[id+3*NSB]*=b*c;
	// 22 flops
#endif
}

__device__ void OrUpdate::setParameter( float param )
{
	orParameter = param;
}

__device__ float OrUpdate::getParameter()
{
	return orParameter;
}
#endif /* ORUPDATE_HXX_ */
