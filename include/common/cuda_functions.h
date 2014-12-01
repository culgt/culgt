/**
 * cuda_functions.h
 *
 *  Created on: Jun 18, 2014
 *      Author: vogt
 */

#ifndef CUDA_FUNCTIONS_H_
#define CUDA_FUNCTIONS_H_

#if !defined(__MATH_FUNCTIONS_H__)

#include <math.h>

inline float rsqrt( float a )
{
	return 1./sqrt( a );
}

inline double rsqrt( double a )
{
	return 1./sqrt( a );
}

#endif


#endif /* CUDA_FUNCTIONS_H_ */
