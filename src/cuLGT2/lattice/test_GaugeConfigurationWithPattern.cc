#include "gmock/gmock.h"
#include <iostream>
#include "GaugeConfiguration.h"
#include "GaugeConfigurationHelper.h"
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
	const LatticeDimension<Ndim> dim;

	static const int mu = 2;
	static const int someIndex = 4;
	static const double someValue;

	LocalLink<SU3Real12<double> > linkWithValue;
	LocalLink<SU3Real12<double> > link;

	typedef StandardPattern<SiteIndex<Ndim,NO_SPLIT>, SU3Real12<double> >  MyPattern;
	typedef GaugeConfiguration<MyPattern> MyGaugeConfigType;
	MyGaugeConfigType* gaugeconfig;

	SiteIndex<Ndim,NO_SPLIT> site;

	GaugeConfigurationWithPattern() : dim(Nt,Nx,Ny,Nz), site(dim){};

	void SetUp()
	{
		site.setLatticeIndex( 3 );
		latticesize = Nt*Nx*Ny*Nz;
		arraysize = Nt*Nx*Ny*Nz*Ndim*12;
		gaugeconfig = new MyGaugeConfigType( dim );
		linkWithValue.set( someIndex, someValue );
	};

	void Teardown()
	{
		gaugeconfig->freeMemory();
		delete gaugeconfig;
		gaugeconfig = NULL;
	};
};

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

TEST_F( GaugeConfigurationWithPattern, SetColdFromGaugeConfigurationHelper )
{
	gaugeconfig->allocateMemoryOnHost();

	GaugeConfigurationHelper<MyPattern>::setCold( gaugeconfig->getHostPointer(), dim );

	ASSERT_FLOAT_EQ( 1.0, gaugeconfig->getElementFromHost(8) );
}
