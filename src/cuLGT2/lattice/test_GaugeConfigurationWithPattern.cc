#include "gmock/gmock.h"
#include <iostream>
#include "GaugeConfiguration.h"
#include "configuration_patterns/StandardPattern.h"
#include "parameterization_types/SU3Real12.h"
#include "../cuLGT1legacy/SiteIndex.hxx"

using namespace culgt;

class GaugeConfigurationWithPattern: public testing::Test {
public:
	int arraysize;
	int latticesize;
	static const int Nc = 3;
	static const int Ndim = 4;
	static const int Nt = 4;
	static const int Nx = 5;
	static const int Ny = 6;
	static const int Nz = 7;
	static const int size[4];

	static const int mu = 2;
	static const int someIndex = 4;
	static const double someValue;

	LocalLink<SU3Real12<double> > linkWithValue;
	LocalLink<SU3Real12<double> > link;

	typedef GaugeConfiguration<StandardPattern<SiteIndex<Ndim,NO_SPLIT>, SU3Real12<double> > > MyGaugeConfigType;
	MyGaugeConfigType* gaugeconfig;

	SiteIndex<Ndim,NO_SPLIT> site;

	GaugeConfigurationWithPattern() : site(size){};

	void SetUp()
	{
		site.setLatticeIndex( 3 );
		latticesize = Nt*Nx*Ny*Nz;
		arraysize = Nt*Nx*Ny*Nz*Ndim*12;
		gaugeconfig = new MyGaugeConfigType( size );
		linkWithValue.set( someIndex, someValue );
	};

	void Teardown()
	{
		gaugeconfig->freeMemory();
		delete gaugeconfig;
		gaugeconfig = NULL;
	};
};

const int GaugeConfigurationWithPattern::size[4] = {Nt,Nx,Ny,Nz};
const double GaugeConfigurationWithPattern::someValue = 13.24;

TEST_F(  GaugeConfigurationWithPattern, SetGetThrowsExceptionIfNoMemoryIsAllocated )
{
	ASSERT_THROW( gaugeconfig->getLinkFromHost( site, mu ), MemoryException );
	ASSERT_THROW( gaugeconfig->setLinkOnHost( site, mu, linkWithValue ), MemoryException );
	ASSERT_THROW( gaugeconfig->getLinkFromDevice( site, mu ), MemoryException );
	ASSERT_THROW( gaugeconfig->setLinkOnDevice( site, mu, linkWithValue ), MemoryException );
}

TEST_F( GaugeConfigurationWithPattern, SetGetLinkOnHost )
{
	gaugeconfig->allocateMemory();

	gaugeconfig->setLinkOnHost( site, mu, linkWithValue );
	link = gaugeconfig->getLinkFromHost( site, mu );

	ASSERT_DOUBLE_EQ( someValue, link.get( someIndex ) );
}

TEST_F( GaugeConfigurationWithPattern, SetLinkOnDevice )
{
	gaugeconfig->allocateMemory();

	gaugeconfig->setLinkOnDevice( site, mu, linkWithValue );
	gaugeconfig->copyToHost();
	link = gaugeconfig->getLinkFromHost( site, mu );

	ASSERT_DOUBLE_EQ( someValue, link.get( someIndex ) );
}

TEST_F( GaugeConfigurationWithPattern, GetLinkFromDevice )
{
	gaugeconfig->allocateMemory();

	gaugeconfig->setLinkOnHost( site, mu, linkWithValue );
	gaugeconfig->copyToDevice();
	link = gaugeconfig->getLinkFromDevice( site, mu );

	ASSERT_DOUBLE_EQ( someValue, link.get( someIndex ) );
}
