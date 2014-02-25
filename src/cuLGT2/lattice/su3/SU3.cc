/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */


/* Random SU3 matrix
 *  0.0956 + 0.5617i   0.2735 + 0.4397i   0.6036 + 0.2069i
 * -0.2839 - 0.4165i  -0.4227 + 0.2138i   0.6040 - 0.3959i
 *  0.6387 - 0.1157i  -0.3653 + 0.6116i  -0.2660 - 0.0218i
 */

#ifndef SU3_H_
#define SU3_H_

#include "../../common/culgt_typedefs.h"

namespace culgt
{

/*
 *
 */
template<typename LinkType> class SU3
{
public:
	LinkType link;

	typename LinkType::PARAMTYPE::TYPE get( lat_group_index_t i )
	{
		return link.get( i );
	}

	void set( lat_group_index_t i, typename LinkType::PARAMTYPE::TYPE val )
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
