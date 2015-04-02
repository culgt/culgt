#include "gmock/gmock.h"

#include "cudatest/cuda_gtest_plugin.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "lattice/GlobalLink.h"
#include "lattice/configuration_patterns/StandardPattern.h"
#include "lattice/configuration_patterns/GPUPattern.h"
#include "lattice/configuration_patterns/GPUPatternParityPriority.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/SU3Real12.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"

using namespace std;

using namespace culgt;
using namespace ::testing;

typedef SiteIndex<4,NO_SPLIT> SiteIndex4;
typedef SiteIndex<4,FULL_SPLIT> SiteIndex4FullSplit;
typedef GlobalLink<StandardPattern<SiteIndex4,SUNRealFull<3,float> > > GlobalLinkStandardSU3Real18Index;
typedef GlobalLink<StandardPattern<SiteIndex4,SU3Real12<float> > > GlobalLinkStandardSU3Real12Index;
typedef GlobalLink<GPUPattern<SiteIndex4,SU3Real12<float> > > GlobalLinkGpuSU3Real12Index;
typedef GlobalLink<GPUPatternParityPriority<SiteIndex4FullSplit,SUNRealFull<3,float> > > GlobalLinkGpuParitySU3Real18Index;

CUDA_TEST(AGlobalLinkWithPatternCu, SetGetValueWithSameLinkAndParameterIndexWorks )
{
	const int size[4]= {2,2,2,2};
	SiteIndex4 s(size, DO_NOT_USE_NEIGHBOURS );
	const int mu = 0;
	const int parameterIndex = 0;
	float U[2*2*2*2*4*18];
	float someValue = 12.421;

	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );

	link.set( parameterIndex, someValue );

	ASSERT_FLOAT_EQ( someValue, link.get( parameterIndex) );
}

CUDA_TEST(AGlobalLinkWithPatternCu, OperatorAssignWithSameTypeDoesNotChangeSite )
{
	const int size[4]= {2,2,2,2};
	SiteIndex4 s(size, DO_NOT_USE_NEIGHBOURS);
	SiteIndex4 s2(size, DO_NOT_USE_NEIGHBOURS);
	const int mu = 0;
	float U[2*2*2*2*4*18];
	int link2Index = 3;

	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );
	s2.setIndex( link2Index );
	GlobalLinkStandardSU3Real18Index link2( U, s2, mu );

	link2 = link;

	ASSERT_EQ( link2Index, link2.site.getIndex() );
}

CUDA_TEST(AGlobalLinkWithPatternCu, OperatorAssignWithSameTypeCopiesData )
{
	const int size[4]= {2,2,2,2};
	SiteIndex4 s(size, DO_NOT_USE_NEIGHBOURS);
	SiteIndex4 s2(size, DO_NOT_USE_NEIGHBOURS);
	const int mu = 0;
	float U[2*2*2*2*4*18];
	float someValue = 12.421;

	s.setIndex( 1 );
	GlobalLinkStandardSU3Real18Index link( U, s, mu );
	s2.setIndex( 3 );
	GlobalLinkStandardSU3Real18Index link2( U, s2, mu );

	link.set( 1, someValue );

	link2.zero();
	link2 = link;

	ASSERT_FLOAT_EQ( someValue, link2.get( 1 ) );
}

CUDA_TEST(AGlobalLinkWithDifferentPatternsCu, OperatorAssignCopiesData )
{
	const int size[4]= {2,2,2,2};
	SiteIndex4 s(size, DO_NOT_USE_NEIGHBOURS);
	SiteIndex4FullSplit s2(size, DO_NOT_USE_NEIGHBOURS);
	const int mu = 0;
	float U[2*2*2*2*4*18];
	float U2[2*2*2*2*4*18];
	float someValue = 12.421;

	s.setIndex( 1 );
	GlobalLinkGpuSU3Real12Index link( U, s, mu );
	s2.setIndex( 3 );
	GlobalLinkGpuParitySU3Real18Index link2( U2, s2, mu );

	link.set( 1, someValue );

	link2.zero();
	link2 = link;

	ASSERT_FLOAT_EQ( someValue, link2.get( 1 ) );
}

