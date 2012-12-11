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
 * Flops: 14
 *
 */

#ifndef MICROUPDATE_HXX_
#define MICROUPDATE_HXX_

#include "../../lattice/datatype/datatypes.h"

class MicroUpdate
{
public:
	__device__ inline MicroUpdate();
	__device__ inline void calculateUpdate( volatile Real (&shA)[4*NSB], short id );
private:
};

__device__ MicroUpdate::MicroUpdate()
{
}

__device__ void MicroUpdate::calculateUpdate( volatile Real (&shA)[4*NSB], short id )
{
#ifdef USE_DP_MICROUPDATE
	double ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	double a0_sq = shA[id]*shA[id];

	double b=2.*shA[id]/(a0_sq+ai_sq);

	shA[id]=(a0_sq-ai_sq)/(a0_sq+ai_sq);
	shA[id+NSB]*=b;
	shA[id+2*NSB]*=b;
	shA[id+3*NSB]*=b;
#else
	Real ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]+shA[id+3*NSB]*shA[id+3*NSB];
	Real a0_sq = shA[id]*shA[id];

	Real b=(Real)2.*shA[id]/(a0_sq+ai_sq);

	shA[id]=(a0_sq-ai_sq)/(a0_sq+ai_sq);
	shA[id+NSB]*=b;
	shA[id+2*NSB]*=b;
	shA[id+3*NSB]*=b;
#endif
}

#endif /* ORUPDATE_HXX_ */
