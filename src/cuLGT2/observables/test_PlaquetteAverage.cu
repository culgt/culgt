

#include "gmock/gmock.h"
#include "../lattice/configuration_patterns/StandardPattern.h"
#include "../lattice/GaugeConfiguration.h"
#include "../lattice/GaugeConfigurationHelper.h"
#include "../cudacommon/DeviceCommunicator.h"
#include "PlaquetteAverage.h"

#include "../cuLGT1legacy/SiteCoord.hxx"
#include "../cuLGT1legacy/SiteIndex.hxx"

using namespace culgt;
using namespace ::testing;

class APlaquetteAverage: public Test
{
public:
	static const int Ndim = 4;
	typedef StandardPattern<SiteIndex<Ndim,NO_SPLIT>, SUNRealFull<3,float> > MyPattern;
	LatticeDimension<Ndim> dim;
	GaugeConfiguration<MyPattern> gc;

//	PlaquetteAverage<MyPattern> plaquetteAverage;

	APlaquetteAverage() : dim( 3,2,2,1 ), gc(dim)
	{
		gc.allocateMemory();
		GaugeConfigurationHelper<MyPattern>::setCold( gc.getHostPointer(), dim );
		gc.copyToDevice();
	}

};

TEST_F( APlaquetteAverage, DeviceCalculatePlaquettesAre1ForColdLattice )
{
	PlaquetteAverage<MyPattern> plaquetteAverage( gc.getDevicePointer(), dim );

	plaquetteAverage.calculatePlaquettes();

	ASSERT_FLOAT_EQ( 1., DeviceCommunicator<float>::getValue( plaquetteAverage.getDevicePointer(), dim.getSize()-1 ) );
}

TEST_F( APlaquetteAverage, ColdLatticeHasAveragePlaquette1 )
{
	PlaquetteAverage<MyPattern> plaquetteAverage( gc.getDevicePointer(), dim );

	ASSERT_FLOAT_EQ( 1.0, plaquetteAverage.getPlaquette() );
}

