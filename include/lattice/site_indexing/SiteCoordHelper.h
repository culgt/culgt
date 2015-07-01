/*
 * SiteCoordHelper.h
 *
 *  Created on: Jul 1, 2015
 *      Author: vogt
 */

#ifndef SITECOORDHELPER_H_
#define SITECOORDHELPER_H_

#include "SiteCoord.h"
#include <vector>
#include "lattice/LatticeDimension.h"

using culgt::SiteCoord;
using std::vector;
using culgt::LatticeDimension;
using culgt::ParityType;

namespace culgt
{

class SiteCoordHelper
{
public:
	template<int Ndim> static SiteCoord<Ndim,NO_SPLIT> makeSite( LatticeDimension<Ndim> dim, vector<int> siteAsVector )
	{
		SiteCoord<Ndim,NO_SPLIT> site( dim );
		for( int i = 0; i < Ndim; i++ )
		{
			site[i] = siteAsVector[i];
		}
		return site;
	}
};

}

#endif /* SITECOORDHELPER_H_ */
