#include "gmock/gmock.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/LocalLink.h"
#include "lattice/GlobalLink.h"
#include "lattice/configuration_patterns/StandardPattern.h"
#include "lattice/configuration_patterns/GPUPattern.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/SU3Real12.h"
#include "lattice/parameterization_types/SU2Real4.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU2_Real4_Real8.h"

using namespace culgt;
using namespace ::testing;

class LocalGlobalLinkWithSU3Real18: public Test
{
public:
	static const int size[4];
	float U[2*2*2*2*4*18];
	SiteIndex<4,NO_SPLIT> site;
	GlobalLink<StandardPattern<SiteIndex<4,NO_SPLIT>,SUNRealFull<3,float> > >* globalLink18;

	LocalLink<SUNRealFull<3,float> > localLink18;

	LocalGlobalLinkWithSU3Real18() : site(size){};

	void SetUp()
	{
		site.setLatticeIndex(1);
		globalLink18 = new GlobalLink<StandardPattern<SiteIndex<4,NO_SPLIT>,SUNRealFull<3,float> > >( U, site, 1 );

		localLink18.zero();
		globalLink18->zero();
	}
};

const int LocalGlobalLinkWithSU3Real18::size[4] = {2,2,2,2};

TEST_F( LocalGlobalLinkWithSU3Real18, OperatorAssignCopiesFromGlobalToLocalForSameParamType )
{
	globalLink18->set( 5, 3.1415 );

	localLink18 = *globalLink18;

	ASSERT_FLOAT_EQ( 3.1415, localLink18.get(5) );
}

TEST_F( LocalGlobalLinkWithSU3Real18, OperatorAssignCopiesFromLocalToGlobalForSameParamType )
{
	localLink18.set( 5, 1.4152 );

	*globalLink18 = localLink18;

	ASSERT_FLOAT_EQ( 1.4152, globalLink18->get(5) );
}


class LocalGlobalLinkWithDifferentParamTypes: public Test
{
public:
	static const int size[4];
	float U[2*2*2*2*4*18];
	SiteIndex<4,NO_SPLIT> site;
	GlobalLink<StandardPattern<SiteIndex<4,NO_SPLIT>,SU3Real12<float> > >* globalLink12;
	LocalLink<SUNRealFull<3,float> > localLink18;

	int someValue;

	LocalGlobalLinkWithDifferentParamTypes() : site( size ), someValue( 1.42 ) {};

	void SetUp()
	{
		site.setLatticeIndex(1);
		globalLink12 = new GlobalLink<StandardPattern<SiteIndex<4,NO_SPLIT>,SU3Real12<float> > >( U, site, 1 );

		localLink18.zero();
		globalLink12->zero();
	}
};
const int LocalGlobalLinkWithDifferentParamTypes::size[4] = {2,2,2,2};

TEST_F( LocalGlobalLinkWithDifferentParamTypes, OperatorAssignCopiesFromGlobal12ToLocal18 )
{
	globalLink12->set( 3, someValue );

	localLink18 = *globalLink12;

	ASSERT_FLOAT_EQ( someValue, localLink18.get( 3 ) );
}

TEST_F( LocalGlobalLinkWithDifferentParamTypes, OperatorAssignCopiesFromGlobal12ToLocal18AndReconstructsThirdLine )
{
	globalLink12->set(0, 1.);
	globalLink12->set(8, 1.);

	localLink18 = *globalLink12;

	ASSERT_FLOAT_EQ( 1., localLink18.get( 16 ) );
}
