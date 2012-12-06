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
