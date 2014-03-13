/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef LOCALLINK_H_
#define LOCALLINK_H_

#include "../cudacommon/cuda_host_device.h"
#include "../common/culgt_typedefs.h"
#include "Link.h"
#include "ParameterizationMediator.h"
#include <cmath>


namespace culgt
{

/**
 * LocalLink allocates memory of size got from ParamType.
 * @author vogt
 */
template<typename ParamType> class LocalLink: Link<ParamType>
{
public:
	/**
	 * Layout of the Link (for example 18 real parameters, 9 complex parameters, or 12 real, ...)
	 */
	typedef ParamType PARAMTYPE;

	CUDA_HOST_DEVICE ~LocalLink()
	{
	};

	/**
	 * Returns the i-th entry of the ParamType specific array
	 * @param i
	 * @return
	 */
	CUDA_HOST_DEVICE inline typename ParamType::TYPE get( lat_group_index_t i ) const
	{
		return store[i];
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	CUDA_HOST_DEVICE inline  void set( lat_group_index_t i, typename ParamType::TYPE val)
	{
		store[i] = val;
	}

	/**
	 * Sets all entries to zero.
	 * @note A function "identity()" (to set the link to the identity) is not possible here, since the identity (for example for SU3) depends on the
	 * 		 Parameterization (ParamType). A delegate should do the job. (Maybe move zero() to the Parameterization class, too.)
	 */
	CUDA_HOST_DEVICE inline  void zero()
	{
		for( lat_group_index_t i = 0; i < ParamType::SIZE; i++ )
		{
			store[i] = (typename ParamType::TYPE) 0.0;
		}
	}

	/**
	 */
	CUDA_HOST_DEVICE inline  void identity()
	{
		ParamType::identity( store );
	}

	/**
	 */
	CUDA_HOST_DEVICE inline typename ParamType::REALTYPE reDet()
	{
		return ParamType::reDet( store );
	}

	CUDA_HOST_DEVICE inline typename ParamType::REALTYPE reTrace()
	{
		return ParamType::reTrace( store );
	}

	/**
	 * Assignment operator for same type just copies all elements.
	 * @param arg
	 * @return
	 */
	CUDA_HOST_DEVICE inline LocalLink<ParamType>& operator=( const LocalLink<ParamType>& arg )
	{
		for( lat_group_index_t i = 0; i < ParamType::SIZE; i++ )
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
	template<typename LinkType> CUDA_HOST_DEVICE inline  LocalLink<ParamType>& operator=( const LinkType arg )
	{
		ParameterizationMediator<ParamType,typename LinkType::PARAMTYPE,LocalLink<ParamType>, LinkType >::assign( *this, arg );
		return *this;
	}
private:
	typename ParamType::TYPE store[ParamType::SIZE];
};


} /* namespace culgt */

#endif /* LOCALLINK_H_ */



