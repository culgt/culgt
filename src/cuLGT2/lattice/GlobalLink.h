/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef GLOBALLINK_H_
#define GLOBALLINK_H_

#include "../common/culgt_typedefs.h"
#include "Link.h"

namespace culgt
{

template<typename ConfigurationPattern> class GlobalLink: Link<typename ConfigurationPattern::PARAMTYPE>
{
public:
	/**
	 * Layout of the Link (for example 18 real parameters, 9 complex parameters, or 12 real, ...)
	 */
	typedef typename ConfigurationPattern::PARAMTYPE PARAMTYPE;
	typedef ConfigurationPattern CONFIGURATIONPATTERN;

	~GlobalLink()
	{
	};

	GlobalLink( typename PARAMTYPE::TYPE* pointerToStore, typename ConfigurationPattern::SITETYPE site, lat_dim_t mu ) : pointerToStore(pointerToStore), site(site), mu(mu) {} ;

	/**
	 * Returns the i-th entry of the ParamType specific array
	 *
	 *  I want something like: store[accesspattern(site,mu,i)]
	 *	where
	 *	- site is an object referring to a lattice site
	 *	- mu is the direction of the link
	 *	- i is the index that belongs to a parameterization
	 *	  (for example with parameterization "float12": i=0 is the real part of the 1/1, i = 1 the imaginary part, i=2 real of 1/2 and so on)
	 * @param i
	 * @return
	 */
	typename PARAMTYPE::TYPE get( lat_group_index_t i ) const
	{
		return pointerToStore[ConfigurationPattern::getIndex(site,mu,i)];
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	void set( lat_group_index_t i, typename PARAMTYPE::TYPE val)
	{
		pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = val;
	}

	typename PARAMTYPE::TYPE* pointerToStore;
	typename ConfigurationPattern::SITETYPE site;
	lat_dim_t mu;
};


} /* namespace culgt */

#endif /* GLOBALLINK_H_ */



