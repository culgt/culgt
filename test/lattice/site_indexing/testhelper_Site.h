/**
 * testhelper_Site.h
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */

#ifndef TESTHELPER_SITE_H_
#define TESTHELPER_SITE_H_

#include <assert.h>

template<typename Site1, typename Site2> bool coordinatesAreEqual( Site1& s1, Site2& s2)
{
	assert(  s1.NDIM == s2.NDIM );
	for( int i = 0; i < s1.NDIM; i++ )
	{
		if( s1[i] != s2[i] ) return false;
	}
	return true;
}


#endif /* TESTHELPER_SITE_H_ */
