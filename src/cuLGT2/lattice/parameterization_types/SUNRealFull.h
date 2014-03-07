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
};

} /* namespace culgt */
#endif /* SUNREALFULL_H_ */
