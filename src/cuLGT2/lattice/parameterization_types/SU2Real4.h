/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU2REAL4_H_
#define SU2REAL4_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

/*
 *
 */
template<typename T> class SU2Real4
{
public:
	static const lat_dim_t NC = 2;
	static const lat_group_index_t SIZE = 4;
	typedef T TYPE;
	typedef T REALTYPE;
};

} /* namespace culgt */
#endif /* SU2REAL4_H_ */
