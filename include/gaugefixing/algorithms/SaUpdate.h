/**
 *
 * This is the Kennedy-Pendleton heatbath update for the SA algorithm
 *
 * In mixed precision (T = float, InternalT = double) the random number are generated as float.
 */

#ifndef SAUPDATE_H_
#define SAUPDATE_H_

#include "../../common/culgt_typedefs.h"
#include "../../util/rng/PhiloxWrapper.h"

namespace culgt
{


template<typename T, typename InternalT=T> class SaUpdate
{
public:
	__device__ inline SaUpdate( float temperature, PhiloxWrapper<T>& rng ) : temperature(temperature), rng(rng){};
	__device__ inline void calculateUpdate( T* shA, const short& id, const int NSB );
	const static int Flops = 86;
private:
	float temperature;
	PhiloxWrapper<T> rng;
};

template<typename T, typename InternalT> __device__ void SaUpdate<T,InternalT>::calculateUpdate( T* shA, const short& id, const int NSB )
{
	InternalT e0,e1,e2,e3, dk, p0;
	InternalT r1,r2,r3,r4;
	InternalT a0,a1,a2,a3;
	InternalT delta, phi, sin_alpha, sin_theta, cos_theta;
	e0=shA[id];
	e1=-shA[id+NSB]; // the minus sign is for the hermitian of the input! be aware of this when reusing this code fragment
	e2=-shA[id+2*NSB]; // "
	e3=-shA[id+3*NSB]; // "
	dk=rsqrt(e0*e0+e1*e1+e2*e2+e3*e3);
	p0=(dk*temperature); // = a*beta

//	16 flop

	do
	{
	  do; while ((r1 = (InternalT)rng.rand2()) < 0.0001);
	  r1 = -::log(r1)*p0;
	  do; while ((r2 = (InternalT)rng.rand2()) < 0.0001);
	  r2 = -::log(r2)*p0;
	  r3 = cospi(2.0*(InternalT)rng.rand2());
	  r3 = r3*r3;
	  delta = r2+r1*r3;
	  r4=(InternalT)rng.rand2();
	} while(r4*r4 > (1.0-0.5*delta));
//	17 flop (if no loops, counting rand() as 1 operation, cospi and sinpi as 2)
	a0=1.0-delta;
	cos_theta=2.0*(InternalT)rng.rand2()-1.0;
	sin_theta=::sqrt(1.0-cos_theta*cos_theta);
	sin_alpha=::sqrt(1-a0*a0);
	phi=2.0*(InternalT)rng.rand2();
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

}

}
#endif
