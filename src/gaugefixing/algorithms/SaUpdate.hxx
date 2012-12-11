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
 * This is the Kennedy-Pendleton heatbath update used in the SA algorithm
 */

#ifndef SAUPDATE_HXX_
#define SAUPDATE_HXX_

#include "../../lattice/datatype/datatypes.h"
#include "../../lattice/rng/PhiloxWrapper.hxx"

class SaUpdate
{
public:
	__device__ inline SaUpdate();
	__device__ inline SaUpdate( float temperature, PhiloxWrapper *rng );
	__device__ inline void calculateUpdate( volatile Real (&shA)[4*NSB], short id );
	__device__ inline void setTemperature( float temperature );
	__device__ inline float getTemperature();
private:
	float temperature;
	PhiloxWrapper *rng;
};

__device__ SaUpdate::SaUpdate()
{
}

__device__ SaUpdate::SaUpdate( float temperature, PhiloxWrapper *rng ) : temperature(temperature), rng(rng)
{
}

__device__ void SaUpdate::calculateUpdate( volatile Real (&shA)[4*NSB], short id )
{
#ifdef USE_DP_SAUPDATE
	// TODO test the DP update in SP code
	double e0,e1,e2,e3, dk, p0;
	double r1,r2,r3,r4;
	double a0,a1,a2,a3;
	double delta, phi, sin_alpha, sin_theta, cos_theta;
	e0=shA[id];
	e1=-shA[id+NSB]; // the minus sign is for the hermitian of the input! be aware of this when reusing this code fragment
	e2=-shA[id+2*NSB]; // "
	e3=-shA[id+3*NSB]; // "
	dk=rsqrt(e0*e0+e1*e1+e2*e2+e3*e3);
	p0=(dk*temperature); // equals a*beta

//	16 flop

	do
	{
	  do; while ((r1 = rng->rand()) < 0.0001);
	  r1 = -log(r1)*p0;
	  do; while ((r2 = rng->rand()) < 0.0001);
	  r2 = -log(r2)*p0;
	  r3 = cospi(2.0*rng->rand());
	  r3 = r3*r3;
	  delta = r2+r1*r3;
	  r4=rng->rand();
	} while(r4*r4 > (1.0-0.5*delta));
//	17 flop (if no loops, counting rand() as 1 operation, cospi/sinpi as 2)
	a0=1.0-delta;
	cos_theta=2.0*rng->rand()-1.0;
	sin_theta=sqrt(1.0-cos_theta*cos_theta);
	sin_alpha=sqrt(1-a0*a0);
	phi=2.0*rng->rand();
	a1=sin_alpha*sin_theta*cospi(phi);
	a2=sin_alpha*sin_theta*sinpi(phi);
	a3=sin_alpha*cos_theta;

	e0 *= dk;
	e1 *= dk;
	e2 *= dk;
	e3 *= dk;
//	25 flop

	shA[id] = a0*e0+a3*e3+a2*e2+e1*a1;
	shA[id+3*NSB] = e0*a3-e3*a0+a1*e2-a2*e1;
	shA[id+2*NSB] = a3*e1-a0*e2+a2*e0-a1*e3;
	shA[id+NSB] = a2*e3+a1*e0-a3*e2-e1*a0;
//	28 flop

//	sum: 86 flop
#else
	Real e0,e1,e2,e3, dk, p0;
	Real r1,r2,r3,r4;
	Real a0,a1,a2,a3;
	Real delta, phi, sin_alpha, sin_theta, cos_theta;
	e0=shA[id];
	e1=-shA[id+NSB]; // the minus sign is for the hermitian of the input! be aware of this when reusing this code fragment
	e2=-shA[id+2*NSB]; // "
	e3=-shA[id+3*NSB]; // "
	dk=rsqrt(e0*e0+e1*e1+e2*e2+e3*e3);
	p0=(dk*temperature); // equals a*beta


	do
	{
	  do; while ((r1 = rng->rand()) < 0.0001);
	  r1 = -log(r1)*p0;
	  do; while ((r2 = rng->rand()) < 0.0001);
	  r2 = -log(r2)*p0;
	  r3 = cospi(2.0*rng->rand());
	  r3 = r3*r3;
	  delta = r2+r1*r3;
	  r4=rng->rand();
	} while(r4*r4 > (1.0-0.5*delta));
	a0=1.0-delta;
	cos_theta=2.0*rng->rand()-1.0;
	sin_theta=sqrt(1.0-cos_theta*cos_theta);
	sin_alpha=sqrt(1-a0*a0);
	phi=2.0*rng->rand();
	a1=sin_alpha*sin_theta*cospi(phi);
	a2=sin_alpha*sin_theta*sinpi(phi);
	a3=sin_alpha*cos_theta;

	e0 *= dk;
	e1 *= dk;
	e2 *= dk;
	e3 *= dk;

	shA[id] = a0*e0+a3*e3+a2*e2+e1*a1;
	shA[id+3*NSB] = e0*a3-e3*a0+a1*e2-a2*e1;
	shA[id+2*NSB] = a3*e1-a0*e2+a2*e0-a1*e3;
	shA[id+NSB] = a2*e3+a1*e0-a3*e2-e1*a0;

#endif
}

__device__ void SaUpdate::setTemperature( float temperature )
{
	this->temperature = temperature;
}

__device__ float SaUpdate::getTemperature()
{
	return temperature;
}
#endif /* ORUPDATE_HXX_ */
