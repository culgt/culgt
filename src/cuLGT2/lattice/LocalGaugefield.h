/**
 * LocalGaugefield.h
 *
 *  Created on: Jun 17, 2014
 *      Author: vogt
 */

#ifndef LOCALGAUGEFIELD_H_
#define LOCALGAUGEFIELD_H_

#include "LinkToGaugefieldConverter.h"

namespace culgt
{

template<typename T, int Nc, gaugefieldtype::GaugefieldType GFType> class LocalGaugefield
{
public:
	template<typename LinkType> LocalGaugefield( LinkType& link )
	{
		LinkToGaugefieldConverter<Nc,GFType>::convert( A, link );
	}

	/**
	 * @param a (1...(Nc*Nc-1)
	 * @return color component
	 */
	T getColorComponent( int a )
	{
		return A[a-1];
	}

	void setColorComponent( int a, T val )
	{
		A[a-1] = val;
	}

private:
	T A[Nc*Nc-1];
};

}

#endif /* LOCALGAUGEFIELD_H_ */
