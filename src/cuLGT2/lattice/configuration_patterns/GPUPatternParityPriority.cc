/**
 * GPUPatternParityPriority.h
 *
 * The pattern is as follows (slowest running index first):
 * parity, mu, paramIndex, siteInParity (where paramIndex for float18 parameterization is a value from 0..17)
 *
 *  Created on: Feb 19, 2014
 *      Author: vogt
 */

#ifndef GPUPATTERNPARITYPRIORITY_H_
#define GPUPATTERNPARITYPRIORITY_H_

#include "../../common/culgt_typedefs.h"
#include "StandardPattern.h"

namespace culgt
{


template<typename Site, typename ParamType> class GPUPatternParityPriority
{
public:
	static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		lat_index_t halfLatSize = site.getSize()/2;
		lat_bool_t parity = site.getIndex() / halfLatSize;
		lat_index_t indexInParity = site.getIndex() % halfLatSize;
		return ((parity*Site::Ndim+mu)*ParamType::SIZE+paramIndex)*halfLatSize+indexInParity;
	}

	/**
	 * For intern-pattern compatibility every pattern should implement a function that converts
	 * its own index to the standard index (which is the index of StandardPattern).
	 * @param index
	 * @param latSize
	 * @return
	 */
	static lat_array_index_t convertToStandardIndex( lat_array_index_t index, lat_index_t latSize )
	{
		lat_index_t halfLatSize = latSize/2;

		lat_index_t siteIndexInParity = index % halfLatSize;

		index /= halfLatSize;
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		index /= ParamType::SIZE;
		lat_dim_t mu = index % Site::Ndim;

		lat_bool_t parity = index / Site::Ndim;

		lat_index_t siteIndex = siteIndexInParity + parity*halfLatSize;

		return StandardPattern<Site, ParamType>::getStandardIndex( siteIndex, mu, paramIndex );
	}
};

}

#endif /* GPUPATTERNPARITYPRIORITY_H_ */


#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"

using namespace culgt;
using namespace ::testing;

TEST( GPUPatternParityPriority, GetIndexReturnsSiteIndexForIndexInLowerHalfOfLattice )
{
	const int paramTypeSize = 10;
	const int latSize = 200;
	const int siteIndex = 10;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int result = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, 0, 0 );

	ASSERT_EQ( mySite.index, result );
}

TEST( GPUPatternParityPriority, GetIndexReturnsSiteIndexForIndexInUpperHalfOfLattice )
{
	const int paramTypeSize = 10;
	const int paramTypeIndex = 1;
	const int latSize = 200;
	const int siteIndex = 110;
	SiteTypeMock<1> mySite(latSize,siteIndex);

	//expect:
	int parity = (siteIndex/(latSize/2)); // = 1
	int siteIndexInParity = (siteIndex%(latSize/2));

	int expect = (parity*paramTypeSize+paramTypeIndex)*(latSize/2)+siteIndexInParity;

	int result = GPUPatternParityPriority<SiteTypeMock<1>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, 0, paramTypeIndex );

	ASSERT_EQ( expect, result );
}


TEST( GPUPatternParityPriority, GetIndexReturnsSiteIndexForIndexInUpperHalfOfLatticeWithMu )
{
	const int paramTypeSize = 10;
	const int paramTypeIndex = 1;
	const int latSize = 200;
	const int siteIndex = 110;
	const int nDim = 4;
	const int myMu = 2;
	SiteTypeMock<nDim> mySite(latSize,siteIndex);

	//expect:
	int parity = (siteIndex/(latSize/2)); // = 1
	int siteIndexInParity = (siteIndex%(latSize/2));

	int expect = ((parity*nDim+myMu)*paramTypeSize+paramTypeIndex)*(latSize/2)+siteIndexInParity;

	int result = GPUPatternParityPriority<SiteTypeMock<nDim>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramTypeIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPatternParityPriority, ConvertToStandardIndexIsCompatibleToStandardPattern )
{
	const int latSize = 200;
	const int siteIndex = 110;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = StandardPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	int gpuPatternIndex = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, myMu, paramIndex );
	int result = GPUPatternParityPriority<SiteTypeMock<>, ParamTypeMock<18> >::convertToStandardIndex( gpuPatternIndex, mySite.getSize() );

	ASSERT_EQ( expect, result );
}
