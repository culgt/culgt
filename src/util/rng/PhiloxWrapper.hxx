/**
 * Wraps the Philox4x32-10 counter based RNG (CBRNG).
 * The benefit of CBRNGs is the absence of large states which are costly to read and write in CUDA code.
 * For a detailed discussion see "Parallel random numbers: as easy as 1, 2, 3" (http://dl.acm.org/citation.cfm?doid=2063405)
 *
 * Any combination of the 4 * 32-bit counter and 2 * 32-bit key gives a uncorrelated random number.
 * As key we use the combination (thread id, seed), where seed can be choosen by the user.
 * The counter is (kernelCounter, globalCounter, *, * ) where
 *  - kernelCounter is a local variable of the kernel and is incremented by each call to generate new random numbers
 *  - globalCounter is a kernel parameter, incrementing has to be done on host side. (TAKE CARE TO DO IT!)
 *    This means that each kernel that uses random numbers has to have an counter parameter which has to be incremented on each kernel call.
 *  - the other two counters are set arbitrarly.
 *
 * Each call to philox4x32 calculates 4 32 bit values.
 *
 * We offer a rand() function that
 *  - returns a Real (float or double)
 *  - takes care to use already calculated numbers.
 * 	- and increments the kernelCounter.
 *
 *
 * TODO: get rid of the preprocessor statements: float/double has to be template argument (in all classes).
 */

#ifndef PHILOXWRAPPER_HXX_
#define PHILOXWRAPPER_HXX_

#include "../datatype/datatypes.h"
#include "../cuda/cuda_host_device.h"
#include "../../external/Random123/philox.h"
#include "../../external/Random123/u01.h"

#include <stdio.h>

class PhiloxWrapper
{
public:
	__device__ inline PhiloxWrapper( int tid, int seed, int globalCounter );
	__device__ inline virtual ~PhiloxWrapper();
	__device__ inline Real rand();

private:
	philox4x32_key_t k;
	philox4x32_ctr_t c;
	union
	{
		philox4x32_ctr_t res;
#ifdef DOUBLEPRECISION
		uint64_t i[2];
#else
		uint32_t i[4];
#endif
	} u;
	short localCounter;
};

PhiloxWrapper::PhiloxWrapper( int tid, int seed, int globalCounter )
{
	k[0] = tid;
	k[1] = seed;
	c[0] = 0; // kernelCounter
	c[1] = globalCounter;
	c[2] = 0x12345678;
	c[3] = 0xabcdef09;
	u.res = philox4x32(c,k);
	localCounter = 0;
}

PhiloxWrapper::~PhiloxWrapper()
{
}

Real PhiloxWrapper::rand()
{
#ifdef DOUBLEPRECISION
	if( localCounter == 0 ) // we have another double available.
	{
		localCounter++;
		return u01_open_open_64_53( u.i[localCounter] );
	}
	else // we have to calculate two new doubles
	{
		localCounter = 0;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_64_53( u.i[localCounter] );
	}
#else
	if( localCounter < 3 ) // we have another float available.
	{
		localCounter++;
		return u01_open_open_32_24( u.i[localCounter] );
	}
	else // we have to calculate 4 new floats
	{
		localCounter = 0;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_32_24( u.i[localCounter] );
	}
#endif
}






#endif /* PHILOXWRAPPER_HXX_ */
