/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef SU3_H_
#define SU3_H_

namespace culgt
{

/*
 *
 */
template<typename LinkType> class SU3
{
public:
	LinkType link;

	typename LinkType::PARAMTYPE::TYPE get( int i )
	{
		return link.get( i );
	}

	void set( int i, typename LinkType::PARAMTYPE::TYPE val )
	{
		link.set( i, val );
	}

	template<typename SU3Type2> SU3<LinkType>& operator=( SU3Type2& val)
	{
		link = val.link;
		return *this;
	}
};

} /* namespace culgt */
#endif /* SU3_H_ */



#include "gmock/gmock.h"
#include "SU3Real12.h"
#include "../LocalLink.h"
#include "SU3Real18.h"

using namespace culgt;
using namespace ::testing;

class ASU3LocalLinkTest: public Test
{
public:
	SU3<LocalLink<SU3Real18<float> > > link;
	SU3<LocalLink<SU3Real18<float> > > link2;
	SU3<LocalLink<SU3Real12<float> > > link3;
};

TEST_F( ASU3LocalLinkTest, CopyingWithSameLinkTypeWorks )
{
	link.set( 2, 1.3 );
	link2 = link;
	ASSERT_FLOAT_EQ( 1.3, link2.get( 2 ) );
}

TEST_F( ASU3LocalLinkTest, CopyingWithDifferentLinkTypeWorks )
{
	link.set( 2, 1.3 );
	link3 = link;
	ASSERT_FLOAT_EQ( 1.3, link3.get( 2 ) );
}
