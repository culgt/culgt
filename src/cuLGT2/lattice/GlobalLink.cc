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

template<typename ConfigurationPattern, typename ParamType> class GlobalLink: Link<ParamType>
{
public:
	/**
	 * Layout of the Link (for example 18 real parameters, 9 complex parameters, or 12 real, ...)
	 */
	typedef ParamType PARAMTYPE;

	~GlobalLink()
	{
	};

	GlobalLink( typename ParamType::TYPE* pointerToStore, typename ConfigurationPattern::SITETYPE site, lat_dim_t mu ) : pointerToStore(pointerToStore), site(site), mu(mu) {} ;

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
	typename ParamType::TYPE get( lat_group_index_t i ) const
	{
		return pointerToStore[ConfigurationPattern::getIndex(site,mu,i)];
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	void set( lat_group_index_t i, typename ParamType::TYPE val)
	{
		pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = val;
	}

private:
	typename ParamType::TYPE* pointerToStore;
	typename ConfigurationPattern::SITETYPE site;
	lat_dim_t mu;
};


} /* namespace culgt */

#endif /* GLOBALLINK_H_ */


#include "gmock/gmock.h"

using namespace culgt;
using namespace ::testing;


class ParameterTypeStub
{
public:
	typedef float TYPE;
};

class SiteStub
{
};

template<int N=0> class AccessPatternStub
{
public:
	typedef SiteStub SITETYPE;

	static int getIndex( SiteStub site, int mu, int patternIndex )
	{
		return N;
	}
};

class AGlobalLink: public Test
{
public:
	SiteStub s;
	static const int mu = 0;
	static const int parameterIndex = 0;
	float U[100];
};

TEST_F(AGlobalLink, SetGetValueWithSameLinkAndParameterIndexWorks )
{
	GlobalLink<AccessPatternStub<>,ParameterTypeStub> link( U, s, mu );
	link.set( parameterIndex, 1.234 );

	ASSERT_FLOAT_EQ( 1.234, link.get( parameterIndex) );
}

TEST_F(AGlobalLink, GetsValueFromTheCorrectPosition )
{
	U[5] = 1.234;
	GlobalLink<AccessPatternStub<5>,ParameterTypeStub> link( U, s, mu );

	ASSERT_FLOAT_EQ( 1.234, link.get(parameterIndex) );
}

TEST_F(AGlobalLink, SetsValueToTheCorrectPosition )
{
	U[5] = -1.;
	GlobalLink<AccessPatternStub<5>,ParameterTypeStub> link( U, s, mu );
	link.set( parameterIndex, 1.234 );

	ASSERT_FLOAT_EQ( 1.234, U[5] );
}



