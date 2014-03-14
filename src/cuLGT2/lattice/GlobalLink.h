/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef GLOBALLINK_H_
#define GLOBALLINK_H_

#include "../cudacommon/cuda_host_device.h"
#include "../common/culgt_typedefs.h"
#include "ParameterizationMediator.h"
#include "Link.h"

namespace culgt
{

template<typename ConfigurationPattern> class GlobalLink//: Link<typename ConfigurationPattern::PARAMTYPE>
{
public:
	/**
	 * Layout of the Link (for example 18 real parameters, 9 complex parameters, or 12 real, ...)
	 */
	typedef typename ConfigurationPattern::PARAMTYPE PARAMTYPE;
	typedef ConfigurationPattern CONFIGURATIONPATTERN;


	CUDA_HOST_DEVICE inline ~GlobalLink()
	{
	};

	CUDA_HOST_DEVICE inline GlobalLink( typename PARAMTYPE::TYPE* pointerToStore, typename ConfigurationPattern::SITETYPE site, lat_dim_t mu ) : pointerToStore(pointerToStore), site(site), mu(mu) {} ;

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
	CUDA_HOST_DEVICE inline typename PARAMTYPE::TYPE get( lat_group_index_t i ) const
	{
		return pointerToStore[ConfigurationPattern::getIndex(site,mu,i)];
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	CUDA_HOST_DEVICE inline void set( lat_group_index_t i, typename PARAMTYPE::TYPE val)
	{
		pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = val;
	}

	CUDA_HOST_DEVICE inline void zero()
	{
		for( lat_group_index_t i = 0; i < PARAMTYPE::SIZE; i++ )
		{
			pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = (typename PARAMTYPE::TYPE) 0.0;
		}
	}

	/**
	 * Assignment operator for same type just copies all elements.
	 * @param arg
	 * @return
	 */
	CUDA_HOST_DEVICE inline GlobalLink<ConfigurationPattern>& operator=( const GlobalLink<ConfigurationPattern>& arg )
	{
		for( lat_group_index_t i = 0; i < ConfigurationPattern::PARAMTYPE::SIZE; i++ )
		{
			set( i, arg.get(i) );
		}
		return *this;
	}

	/**
	 * Assignment operator needs to call Mediator for different ParamTypes and/or LinkTypes
	 * @param arg
	 * @return
	 */
	template<typename LinkType> inline CUDA_HOST_DEVICE GlobalLink<ConfigurationPattern>& operator=( const LinkType arg )
	{
		ParameterizationMediator<typename ConfigurationPattern::PARAMTYPE,typename LinkType::PARAMTYPE,GlobalLink<ConfigurationPattern>, LinkType >::assign( *this, arg );
		return *this;
	}

	typename PARAMTYPE::TYPE* pointerToStore;
	typename ConfigurationPattern::SITETYPE site;
	lat_dim_t mu;
};


} /* namespace culgt */

#endif /* GLOBALLINK_H_ */



