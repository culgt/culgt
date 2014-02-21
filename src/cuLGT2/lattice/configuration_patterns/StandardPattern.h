/**
 * StandardPattern.h
 *
 * The pattern is as follows (slowest running index first):
 * site, mu, paramIndex (where paramIndex for float18 parameterization is a value from 0..17)
 *
 *  Created on: Feb 19, 2014
 *      Author: vogt
 */

#ifndef STANDARDPATTERN_H_
#define STANDARDPATTERN_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

template<typename Site, typename ParamType> class StandardPattern
{
public:
	typedef Site SITETYPE;
	typedef ParamType PARAMTYPE;

	static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		return (site.getIndex()*Site::Ndim+mu)*ParamType::SIZE+paramIndex;
	}

	/**
	 * Use in other patterns to compute the StandardIndex from siteIndex, mu, paramIndex
	 * @param siteIndex
	 * @param mu
	 * @param paramIndex
	 * @param latSize
	 * @return
	 */
	static lat_array_index_t getStandardIndex( lat_index_t siteIndex, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		return (siteIndex*Site::Ndim+mu)*ParamType::SIZE+paramIndex;
	}
};

}

#endif /* STANDARDPATTERN_H_ */


