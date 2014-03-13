/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SUNREALFULL_H_
#define SUNREALFULL_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

/*
 *
 */
template<int Nc, typename T> class SUNRealFull
{
public:
	static const lat_dim_t NC = Nc;
	static const lat_group_index_t SIZE = 2*(Nc*Nc);
	typedef T TYPE;
	typedef T REALTYPE;

	static void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
		for( int i = 0; i < Nc; i++ )
		{
			store[2*i*(1+Nc)] = 1.0;
		}
	}

	static REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		REALTYPE result = 0;
		for( int i = 0; i < Nc; i++ )
		{
			result += store[2*i*(1+Nc)];
		}
		return result;
	}
};

template<typename T> class SUNRealFull<3,T>
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = 2*(3*3);
	typedef T TYPE;
	typedef T REALTYPE;

	static void inline identity( TYPE store[SIZE] )
	{
		for( int i = 0; i < SIZE; i++ )
		{
			store[i] = 0;
		}
		for( int i = 0; i < 3; i++ )
		{
			store[2*i*(1+3)] = 1.0;
		}
	}

	static inline REALTYPE reDet( TYPE store[SIZE] )
	{
		return store[0]*store[8]*store[16] - store[0]*store[10]*store[14] - store[2]*store[6]*store[16] + store[2]*store[10]*store[12] + store[4]*store[6]*store[14] - store[4]*store[8]*store[12] - store[0]*store[9]*store[17] + store[0]*store[11]*store[15] - store[1]*store[8]*store[17] - store[1]*store[9]*store[16] + store[1]*store[10]*store[15] + store[1]*store[11]*store[14] + store[2]*store[7]*store[17] - store[2]*store[11]*store[13] + store[3]*store[6]*store[17] + store[3]*store[7]*store[16] - store[3]*store[10]*store[13] - store[3]*store[11]*store[12] - store[4]*store[7]*store[15] + store[4]*store[9]*store[13] - store[5]*store[6]*store[15] - store[5]*store[7]*store[14] + store[5]*store[8]*store[13] + store[5]*store[9]*store[12];
	}

	static REALTYPE inline reTrace( TYPE store[SIZE] )
	{
		return store[0]+store[8]+store[16];
	}
};


} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
