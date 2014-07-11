/**
 * mathfunction_overload.h
 *
 *  Created on: Jul 7, 2014
 *      Author: vogt
 */

#ifndef MATHFUNCTION_OVERLOAD_H_
#define MATHFUNCTION_OVERLOAD_H_
#include "../cudacommon/cuda_host_device.h"

namespace culgt
{
	inline double norm_squared( double& a )
	{
		return a*a;
	}
	inline float norm_squared( float& a )
	{
		return a*a;
	}

	inline double conj( double& a )
	{
		return a;
	}
	inline float conj( float& a )
	{
		return a;
	}

	template<typename T> CUDA_HOST_DEVICE T sign( T a )
	{
		if( a > 0 ) return 1.;
		else if( a < 0 ) return -1.;
		else return 0.;
	}


//	float fabs( float a )
//	{
//		return std::fabs( a );
//	}
//	double fabs( double a )
//	{
//		return std::fabs( a );
//	}

//	double cos( double a )
//	{
//		return std::cos( a );
//	}
//	float cos( float a )
//	{
//		return std::cos( a );
//	}
//
//	double acos( double a )
//	{
//		return std::acos( a );
//	}
//	float acos( float a )
//	{
//		return std::acos( a );
//	}
//
//	double sqrt( double a )
//	{
//		return std::sqrt( a );
//	}
//	float sqrt( float a )
//	{
//		return std::sqrt( a );
//	}


}


#endif /* MATHFUNCTION_OVERLOAD_H_ */
