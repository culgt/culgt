
#include "gmock/gmock.h"
#include "lattice/GaugeConfiguration.h"
#include "lattice/parameterization_types/SU2Vector4.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "lattice/configuration_patterns/GPUPatternParityPriority.h"
#include "lattice/filetypes/LinkFileMDP.h"
#include "observables/PlaquetteAverage.h"
using namespace culgt;
using namespace ::testing;

template<int N> class LinkFileMDPSU2Template: public Test
{
public:
	typedef double REAL;
	typedef SU2Vector4<REAL> PARAMTYPE;
	typedef SiteIndex<4,FULL_SPLIT> SITE;
	typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
	typedef LocalLink<SUNRealFull<2,REAL> > LOCALLINK;

	LatticeDimension<4> dim;
	GaugeConfiguration<PATTERNTYPE> config;
	GaugeConfiguration<PATTERNTYPE> config2;
	PlaquetteAverage<PATTERNTYPE,LOCALLINK>* plaquetteCalculator;

	const double plaquetteValue = 0.58416374374822344;

	void SetUp()
	{
		config.allocateMemory();
		DeviceMemoryManager::clear( config.getDevicePointer() );
		config2.allocateMemory();
		plaquetteCalculator = new PlaquetteAverage<PATTERNTYPE,LOCALLINK>( config.getDevicePointer(), dim );
	}
	void TearDown()
	{
		config.freeMemory();
		config2.freeMemory();
		delete plaquetteCalculator;
	}

	REAL calcPlaquetteOnConfig()
	{
		return plaquetteCalculator->getPlaquette();
	}

	LinkFileMDPSU2Template(): dim(N,N,N,N), config( dim ), config2(dim)
	{
	}
};

typedef LinkFileMDPSU2Template<4> LinkFileMDPSU2;

TEST_F( LinkFileMDPSU2, CheckPlaquette )
{
	LinkFileMDP<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "sample_4x4x4x4_su2.orig.mdp" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_DOUBLE_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileMDPSU2, CheckPlaquetteAfterReadWriteRead )
{
	LinkFileMDP<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "sample_4x4x4x4_su2.orig.mdp" );
	config2.loadFile( linkfileIn );

	LinkFileMDP<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "temp.mdp" );
	config2.saveFile( linkfileOut );

	LinkFileMDP<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "temp.mdp" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_DOUBLE_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}




