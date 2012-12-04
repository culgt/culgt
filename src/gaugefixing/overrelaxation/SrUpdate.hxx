/*
 * SrUpdate.hxx: stochastic relaxation
 *
 *  Created on: June 20, 2012
 *      Author: schroeck
 */

#ifndef SRUPDATE_HXX_
#define SRUPDATE_HXX_

#include "../../../util/rng/PhiloxWrapper.hxx"
#include "../../../util/datatype/datatypes.h"

class SrUpdate
{
public:
	__device__ inline SrUpdate();
	__device__ inline SrUpdate( float param, PhiloxWrapper *rng );
	__device__ inline void calculateUpdate( volatile Real (&shA)[128], short id );
	__device__ inline void setParameter( float param );
	__device__ inline float getParameter();
private:
	float srParameter;
	PhiloxWrapper *rng;
};

__device__ SrUpdate::SrUpdate()
{
}

__device__ SrUpdate::SrUpdate( float param, PhiloxWrapper *rng  ) : srParameter(param), rng(rng)
{
}

__device__ void SrUpdate::calculateUpdate( volatile Real (&shA)[128], short id )
{
//stochastic overrelaxation
#ifdef USE_DP_SRUPDATE
	double rand = rng->rand();
	double a0,a1,a2,a3,c;
	a0 = shA[id];
	a1 = shA[id+32];
	a2 = shA[id+64];
	a3 = shA[id+96];

	shA[id]    = (rand>=srParameter)*a0 + (rand<srParameter)*(a0*a0-a1*a1-a2*a2-a3*a3); // 12 flop
	shA[id+32] = (rand>=srParameter)*a1 + (rand<srParameter)*(2.0*a0*a1); // 7 flop
	shA[id+64] = (rand>=srParameter)*a2 + (rand<srParameter)*(2.0*a0*a2);
	shA[id+96] = (rand>=srParameter)*a3 + (rand<srParameter)*(2.0*a0*a3);

	c=rsqrt(shA[id]*shA[id]+shA[id+32]*shA[id+32]+shA[id+64]*shA[id+64]+shA[id+96]*shA[id+96]); // 8 flop
	shA[id]    *= c;
	shA[id+32] *= c;
	shA[id+64] *= c;
	shA[id+96] *= c;
	// 4 flop

	// sum: 31 flop
#else
	Real rand = rng->rand();
	Real a0,a1,a2,a3,c;
	a0 = shA[id];
	a1 = shA[id+32];
	a2 = shA[id+64];
	a3 = shA[id+96];
	
	shA[id]    = (rand>=srParameter)*a0 + (rand<srParameter)*(a0*a0-a1*a1-a2*a2-a3*a3);
	shA[id+32] = (rand>=srParameter)*a1 + (rand<srParameter)*(2.0*a0*a1);
	shA[id+64] = (rand>=srParameter)*a2 + (rand<srParameter)*(2.0*a0*a2);
	shA[id+96] = (rand>=srParameter)*a3 + (rand<srParameter)*(2.0*a0*a3);

	c=rsqrt(shA[id]*shA[id]+shA[id+32]*shA[id+32]+shA[id+64]*shA[id+64]+shA[id+96]*shA[id+96]);
	shA[id]    *= c;
	shA[id+32] *= c;
	shA[id+64] *= c;
	shA[id+96] *= c;
#endif
}

__device__ void SrUpdate::setParameter( float param )
{
	srParameter = param;
}

__device__ float SrUpdate::getParameter()
{
	return srParameter;
}
#endif /* SRUPDATE_HXX_ */
