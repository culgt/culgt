/**
 * GPUPattern.h
 *
 * The pattern is as follows (slowest running index first):
 * mu, paramIndex, site (where paramIndex for float18 parameterization is a value from 0..17)
 *
 * TODO Test difference to paramIndex,mu,site
 *
 *  Created on: Feb 19, 2014
 *      Author: vogt
 */

#ifndef GPUPATTERN_H_
#define GPUPATTERN_H_

#include "../../common/culgt_typedefs.h"

#include "StandardPattern.h"
namespace culgt
{


template<typename Site, typename ParamType> class GPUPattern
{
public:
	static lat_array_index_t getIndex( const Site& site, lat_dim_t mu, lat_group_index_t paramIndex )
	{
		return (mu*ParamType::SIZE+paramIndex)*site.getSize()+site.getIndex();
	}

	/**
	 * For inter-pattern compatibility every pattern should implement a function that converts
	 * its own index to the standard index (which is the index of StandardPattern).
	 * @param index
	 * @param latSize
	 * @return
	 */
	static lat_array_index_t convertToStandardIndex( lat_array_index_t index, lat_index_t latSize )
	{
		lat_index_t siteIndex = index % latSize;

		index /= latSize;
		lat_group_index_t paramIndex = index % ParamType::SIZE;

		lat_dim_t mu = index / ParamType::SIZE;

		return StandardPattern<Site, ParamType>::getStandardIndex( siteIndex, mu, paramIndex );
	}
};

}

#endif /* GPUPATTERN_H_ */


#include "gmock/gmock.h"
#include "testhelper_PatternMocks.h"

using namespace culgt;
using namespace ::testing;

TEST( GPUPattern, GetIndexReturnsSiteIndexForParamIndexZeroAndMuZero )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<0> >::getIndex( mySite, 0, 0 );

	ASSERT_EQ( mySite.index, result );
}

TEST( GPUPattern, GetIndexReturnsCorrectIndexWithMuZero )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramIndex = 4;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = paramIndex*latSize+siteIndex;

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, 0, paramIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPattern, GetIndexReturnsCorrectIndex )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = (myMu*paramTypeSize+paramIndex)*latSize+siteIndex;

	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	ASSERT_EQ( expect, result );
}

TEST( GPUPattern, ConvertToStandardIndexIsCompatibleToStandardPattern )
{
	const int latSize = 1000;
	const int siteIndex = 100;
	const int paramTypeSize = 18;
	const int paramIndex = 4;
	const int myMu = 3;
	const SiteTypeMock<> mySite(latSize,siteIndex);

	int expect = StandardPattern<SiteTypeMock<>, ParamTypeMock<paramTypeSize> >::getIndex( mySite, myMu, paramIndex );

	int gpuPatternIndex = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::getIndex( mySite, myMu, paramIndex );
	int result = GPUPattern<SiteTypeMock<>, ParamTypeMock<18> >::convertToStandardIndex( gpuPatternIndex, mySite.getSize() );

	ASSERT_EQ( expect, result );
}
