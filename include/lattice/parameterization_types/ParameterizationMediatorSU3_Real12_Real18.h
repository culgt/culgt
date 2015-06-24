/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_
#define PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_

#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"
#include "ParameterizationMediator.h"
#include "SU3Real12.h"
#include "SUNRealFull.h"

namespace culgt
{

template<typename LinkType1, typename LinkType2, typename T> class ParameterizationMediator<SU3Real12<T>,SUNRealFull<3,T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 12; i++ )
		{
			l1.set(i, l2.get(i));
		}
	}
};

//template<typename LinkType1, typename LinkType2> class ParameterizationMediator<SU3Real12<float>,SUNRealFull<3,float>, LinkType1, LinkType2>
//{
//public:
//	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
//	{
//		for( lat_group_index_t i = 0; i < 12; i++ )
//		{
//			l1.set(i, l2.get(i));
//		}
//	}
//};

/**
 * Copy from 12 to 18 parameterization needs reconstruction of "third row"
 */
template<typename LinkType1, typename LinkType2, typename T>class ParameterizationMediator<SUNRealFull<3,T>,SU3Real12<T>, LinkType1, LinkType2>
{
public:
	CUDA_HOST_DEVICE inline static void assign( LinkType1& l1, const LinkType2& l2 )
	{
		for( lat_group_index_t i = 0; i < 12; i++ )
		{
			l1.set(i, l2.get(i));
		}

		T a12 = l1.get(2)*l1.get(10)-l1.get(3)*l1.get(11) - l1.get(8)*l1.get(4)+l1.get(9)*l1.get(5);
		T a13 = -l1.get(2)*l1.get(11)-l1.get(3)*l1.get(10) + l1.get(8)*l1.get(5)+l1.get(9)*l1.get(4);
		T a14 = l1.get(4)*l1.get(6)-l1.get(5)*l1.get(7) - l1.get(10)*l1.get(0)+l1.get(11)*l1.get(1) ;
		T a15 = -l1.get(4)*l1.get(7)-l1.get(5)*l1.get(6) + l1.get(10)*l1.get(1)+l1.get(11)*l1.get(0);
		T a16 = l1.get(0)*l1.get(8)-l1.get(1)*l1.get(9) - l1.get(2)*l1.get(6)+l1.get(3)*l1.get(7);
		T a17 = -l1.get(0)*l1.get(9)-l1.get(1)*l1.get(8) + l1.get(2)*l1.get(7)+l1.get(3)*l1.get(6);

		T norm = a12*a12+a13*a13+a14*a14+a15*a15+a16*a16+a17*a17;
		norm = rsqrt( norm );

		l1.set(12, a12*norm );
		l1.set(13, a13*norm );
		l1.set(14, a14*norm );
		l1.set(15, a15*norm );
		l1.set(16, a16*norm );
		l1.set(17, a17*norm );
	}
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_REAL12_REAL18_H_ */



