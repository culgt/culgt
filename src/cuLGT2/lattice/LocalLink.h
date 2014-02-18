/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef LOCALLINK_H_
#define LOCALLINK_H_

#include "Link.h"
#include "ParameterizationMediator.h"

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

	~LocalLink()
	{
	};

	/**
	 * Returns the i-th entry of the ParamType specific array
	 * @param i
	 * @return
	 */
	typename ParamType::TYPE get( int i ) const
	{
		return store[i];
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	void set( int i, typename ParamType::TYPE val)
	{
		store[i] = val;
	}

	/**
	 * Sets all entries to zero.
	 * @note A function "identity()" is not possible here, since the identity (for example for SU3) depends on the
	 * 		 Parameterization (ParamType).
	 */
	void zero()
	{
		for( int i = 0; i < ParamType::SIZE; i++ )
		{
			store[i] = (typename ParamType::TYPE) 0.0;
		}
	}

	/**
	 * Assignment operator for same type just copies all elements.
	 * @param arg
	 * @return
	 */
	LocalLink<ParamType>& operator=( const LocalLink<ParamType>& arg )
	{
		for( int i = 0; i < ParamType::SIZE; i++ )
		{
			set( i, arg.get(i) );
		}
		return *this;
	}

	/**
	 * Assignment operator for different types needs to delegate the work to a mediator that
	 * has precise knowledge of both parameterizations.
	 * @param arg
	 * @return
	 */
	template<typename ParamType2> LocalLink<ParamType>& operator=( const LocalLink<ParamType2>& arg )
	{
		ParameterizationMediator<ParamType,ParamType2,LocalLink<ParamType>, LocalLink<ParamType2> >::assign( *this, arg );
		return *this;
	}
private:
	typename ParamType::TYPE store[ParamType::SIZE];
};


} /* namespace culgt */

#endif /* LOCALLINK_H_ */



