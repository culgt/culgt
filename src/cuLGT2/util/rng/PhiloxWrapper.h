/**
 *
 * Wraps the Philox4x32-10 counter based RNG (CBRNG).
 * The benefit of CBRNGs is the absence of large states which are costly to read and write in CUDA code.
 * For a detailed discussion see "Parallel random numbers: as easy as 1, 2, 3" (http://dl.acm.org/citation.cfm?doid=2063405)
 *
 * Any combination of the 4 * 32-bit counter and 2 * 32-bit key gives a uncorrelated random number.
 * As key we use the combination (thread id, seed), where seed can be chosen by the user.
 * The counter is (kernelCounter, globalCounter, *, * ) where
 *  - kernelCounter is a local variable of the kernel and is incremented by each call to generate new random numbers
 *  - globalCounter is a kernel parameter, incrementing has to be done on host side. (TAKE CARE TO DO IT!)
 *    This means that each kernel that uses random numbers has to have an counter parameter which has to be incremented on each kernel call.
 *  - the other two counters are set arbitrarily.
 *
 * Each call to philox4x32 calculates 4 32 bit values.
 *
 * We offer a rand() function that
 *  - returns a Real (float or double)
 *  - takes care to use already calculated numbers.
 * 	- and increments the kernelCounter.
 *
 * The static getNextCounter() function gives a global (runtime-wide) counter which can be given to the constructor.
 * I don't want to do this implicitly to avoid mixing of host and device variables.
 *
 *
 * TODO: get rid of the preprocessor statements: float/double has to be template argument (in all classes).
 * 		Because now there is no code possible that uses both single and double precision random numbers.
 * TODO: A lot of stack-frame is used (at least in in CUDA5.0) when Philox is involved (~ 200 bytes for the single precision SA kernel)
 */

#ifndef PHILOXWRAPPER_H_
#define PHILOXWRAPPER_H_

#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"
#include "../../external/Random123/philox.h"
#include "../../external/Random123/u01.h"

#include <cuda.h>

#include <stdio.h>

namespace culgt
{


template<class T> struct PhiloxWrapperInfo
{
};

template<> struct PhiloxWrapperInfo<float>
{
    static const int uintSize = 4;
    typedef uint32_t uintType;
    typedef float2 T2;
};
template<> struct PhiloxWrapperInfo<double>
{
    static const int uintSize = 2;
    typedef uint64_t uintType;
    typedef double2 T2;
};



template<typename T> class PhiloxWrapper
{
public:
	__device__ inline PhiloxWrapper( int tid, int seed, unsigned int globalCounter );
	__device__ inline virtual ~PhiloxWrapper();
	__device__ inline T rand();
	__device__ inline T rand2();
	__device__ inline T randExp( T lambda);
	__device__ inline typename PhiloxWrapperInfo<T>::T2 randNormal();

	static __host__ unsigned int getNextCounter();
	static __host__ unsigned int getCurrentCounter();

private:
	philox4x32_key_t k;
	philox4x32_ctr_t c;
	union
	{
		philox4x32_ctr_t res;
		typename PhiloxWrapperInfo<T>::uintType i[PhiloxWrapperInfo<T>::uintSize];
	} u;
	short localCounter;
	static unsigned int globalCounter;
};


template<typename T> unsigned int PhiloxWrapper<T>::globalCounter = 0;

template<typename T> __device__ PhiloxWrapper<T>::PhiloxWrapper( int tid, int seed, unsigned int globalCounter )
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

template<typename T> __device__ PhiloxWrapper<T>::~PhiloxWrapper()
{
}

template<> __device__ double PhiloxWrapper<double>::rand()
{
	if( localCounter == 0 )// we have to calculate two new doubles
	{
		localCounter = 1;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_64_53( u.i[localCounter] );
	}
	else  // we have another double available.
	{
		localCounter--;
		return u01_open_open_64_53( u.i[localCounter] );
	}
}

template<> __device__ float PhiloxWrapper<float>::rand()
{
	if( localCounter == 0 ) // we have to calculate 4 new floats
	{
		localCounter = 3;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_32_24( u.i[localCounter] );
	}
	else // we have another float available.
	{
		localCounter--;
		return u01_open_open_32_24( u.i[localCounter] );
	}
}

/**
 * Tested in Cuda5.5: This method uses no stack frame (compared to 52 bytes of rand()) but 30 registers (rand() can be forced to use 23 registers with maxrregcount otherwise it will use 31 but with better performance).
 * Still it is slower. Maybe for some applications this method will be better. (Maybe in case of register pressure?)
 *
 * I also tried a variant which generates a new set after every 2 random numbers (instead of 4 for float). As expect, this is also slightly slower. The idea was inspired by the (comparably) very fast double precision rngs!
 * TODO: check the philox2x32
 */
template<> __device__ float PhiloxWrapper<float>::rand2()
{
	switch( localCounter )
	{
	case 0:
		localCounter = 1;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_32_24( u.i[localCounter] );
	case 1:
		localCounter = 0;
		return u01_open_open_32_24( u.i[localCounter] );
	case 2:
		localCounter = 1;
		return u01_open_open_32_24( u.i[localCounter] );
	default:
		localCounter = 2;
		return u01_open_open_32_24( u.i[localCounter] );
	}
}

template<> __device__ double PhiloxWrapper<double>::rand2()
{
	switch( localCounter )
	{
	case 0:
		localCounter = 1;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return u01_open_open_64_53( u.i[localCounter] );
	default:
		localCounter = 0;
		return u01_open_open_64_53( u.i[localCounter] );
	}
}

template<> __device__ double PhiloxWrapper<double>::randExp( double lambda)
{
	if( localCounter == 0 )// we have to calculate two new doubles
	{
		localCounter = 1;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return -log(u01_open_closed_64_53( u.i[localCounter] ))/lambda;
	}
	else  // we have another double available.
	{
		localCounter--;
		return -log(u01_open_closed_64_53( u.i[localCounter] ))/lambda;
	}
}

template<> __device__ float PhiloxWrapper<float>::randExp( float lambda )
{
	if( localCounter == 0 ) // we have to calculate 4 new floats
	{
		localCounter = 3;
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);
		return -log(u01_open_closed_32_24( u.i[localCounter] ))/lambda;
	}
	else // we have another float available.
	{
		localCounter--;
		return -log(u01_open_closed_32_24( u.i[localCounter] ))/lambda;
	}
}

/**
 * Returns two standard normal distributed random numbers using Box-Muller:
 *
 *
 *
 * @return
 */
template<> __device__ double2 PhiloxWrapper<double>::randNormal()
{
	// we calculate 2 new random numbers
	localCounter = 0; // ensure that the next call uses new random numbers
	c[0]++;
	u.res = philox4x32(c,k);

	double2 normal;

	double temp = sqrt( -2.*log(u01_open_open_64_53(u.i[0])));
	double temp2 = u01_open_open_64_53( u.i[1] );
	normal.x = temp*cospi(2.*temp2);
	normal.y = temp*sinpi(2.*temp2);

	return normal;
}

template<> __device__ float2 PhiloxWrapper<float>::randNormal()
{
	// if only one number is still available: discard and generate 4 new

	float2 normal;
	if( localCounter <= 1 ) // we have to calculate new floats
	{
		localCounter = 2; // we use already 2 numbers
		c[0]++; // inkrement the kernel counter
		u.res = philox4x32(c,k);

		double temp = sqrt( -2.*log(u01_open_open_32_24(u.i[3])));
		double temp2 = u01_open_open_32_24(u.i[2]);
		normal.x = temp*cospi(2.*temp2);
		normal.y = temp*sinpi(2.*temp2);

	}
	else // we have another float available.
	{
		localCounter--;
		double temp = sqrt( -2.*log(u01_open_open_32_24(u.i[localCounter])));
		localCounter--;
		double temp2 = u01_open_open_32_24(u.i[localCounter]);
		normal.x = temp*cospi(2.*temp2);
		normal.y = temp*sinpi(2.*temp2);
	}

	return normal;
}

/**
 * Use this to get a global (static) counter to initialize the kernel-specific PhiloxWrapper.
 */
template<typename T> __host__ unsigned int PhiloxWrapper<T>::getNextCounter()
{
	return globalCounter++;
}

template<typename T> __host__ unsigned int PhiloxWrapper<T>::getCurrentCounter()
{
	return globalCounter;
}

}

#endif /* PHILOXWRAPPER_HXX_ */
