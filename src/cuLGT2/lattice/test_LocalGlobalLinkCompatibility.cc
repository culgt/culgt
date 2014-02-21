#include "gmock/gmock.h"
#include "LocalLink.h"
#include "GlobalLink.h"
#include "configuration_patterns/StandardPattern.h"
#include "../cuLGT1legacy/SiteIndex.hxx"
#include "su3/SU3Real12.h"
#include "su3/ParameterizationMediatorSU3_Real12_Real18.h"

using namespace culgt;
using namespace ::testing;

TEST( LocalLinkGlobalLinkWithSU3Real18, OperatorAssignCopiesFromGlobalToLocalForSameParamType )
{
	LocalLink<SU3Real18<float> > localLink18;
	localLink18.zero();

	int size[4] = {2,2,2,2};
	float U[2*2*2*2*4*18];
	SiteIndex<4,NO_SPLIT> site( size );
	site.setLatticeIndex(1);
	GlobalLink<StandardPattern<SiteIndex<4,NO_SPLIT>,SU3Real18<float> > > globalLink18( U, site, 1 );

	localLink18 = globalLink18;
}

