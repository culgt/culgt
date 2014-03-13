/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU3REAL12_H_
#define SU3REAL12_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

/*
 *
 */
template<typename T> class SU3Real12
{
public:
	static const lat_dim_t NC = 3;
	static const lat_group_index_t SIZE = 12;
	typedef T TYPE;
	typedef T REALTYPE;
};

} /* namespace culgt */
#endif /* SU3REAL12_H_ */
