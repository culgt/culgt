/*
 * Overrelaxation.hxx
 *
 *  Created on: May 10, 2012
 *      Author: vogt
 *
 *      14 Flops
 */

#ifndef MICROUPDATE_HXX_
#define MICROUPDATE_HXX_

#include "../../../util/datatype/datatypes.h"

class MicroUpdate
{
public:
	__device__ inline MicroUpdate();
	__device__ inline void calculateUpdate( volatile Real (&shA)[NSB4], short id );
private:
};

__device__ MicroUpdate::MicroUpdate()
{
}

__device__ void MicroUpdate::calculateUpdate( volatile Real (&shA)[NSB4], short id )
{
#ifdef USE_DP_MICROUPDATE
	double ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+NSB2]*shA[id+NSB2]+shA[id+NSB3]*shA[id+NSB3];
	double a0_sq = shA[id]*shA[id];

	double b=2.*shA[id]/(a0_sq+ai_sq);

	shA[id]=(a0_sq-ai_sq)/(a0_sq+ai_sq);
	shA[id+NSB]*=b;
	shA[id+NSB2]*=b;
	shA[id+NSB3]*=b;
#else
	Real ai_sq = shA[id+NSB]*shA[id+NSB]+shA[id+NSB2]*shA[id+NSB2]+shA[id+NSB3]*shA[id+NSB3];
	Real a0_sq = shA[id]*shA[id];

	Real b=(Real)2.*shA[id]/(a0_sq+ai_sq);

	shA[id]=(a0_sq-ai_sq)/(a0_sq+ai_sq);
	shA[id+NSB]*=b;
	shA[id+NSB2]*=b;
	shA[id+NSB3]*=b;
#endif
}

#endif /* ORUPDATE_HXX_ */
