/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef LINK_H_
#define LINK_H_

namespace culgt
{

/*
 *
 */
template<typename ParamType> class Link
{
public:
	virtual typename ParamType::TYPE get( lat_group_index_t i ) const = 0;
	virtual void set( lat_group_index_t i, typename ParamType::TYPE val ) = 0;
	virtual ~Link(){};
};

} /* namespace culgt */
#endif /* LINK_H_ */
