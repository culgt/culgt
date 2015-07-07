/**
 * FIXME: should not need cuda (needs for calculating plaquette)
 */

#include "gmock/gmock.h"
#include "lattice/GaugeConfiguration.h"
#include "lattice/parameterization_types/SU3Vector4.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "lattice/configuration_patterns/GPUPatternParityPriority.h"
#include "lattice/filetypes/filetype_config.h"
#include "lattice/filetypes/LinkFileVogt.h"
#include "lattice/filetypes/LinkFileNERSC.h"
#include "lattice/filetypes/LinkFileMDP.h"
#include "lattice/filetypes/LinkFileHirep.h"
#include "observables/PlaquetteAverage.h"

#ifdef CULGT_HAVE_LINKFILE_ILDG
#include "lattice/filetypes/LinkFileILDG.h"
#endif

using namespace culgt;
using namespace ::testing;

template<int N> class LinkFileCompatibilitySU3Template: public Test
{
public:
	typedef float REAL;
	typedef SU3Vector4<REAL> PARAMTYPE;
	typedef SiteIndex<4,FULL_SPLIT> SITE;
	typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
	typedef LocalLink<SUNRealFull<3,REAL> > LOCALLINK;

	LatticeDimension<4> dim;
	GaugeConfiguration<PATTERNTYPE> config;
	GaugeConfiguration<PATTERNTYPE> config2;
	PlaquetteAverage<PATTERNTYPE,LOCALLINK>* plaquetteCalculator;

	const double plaquetteValue = 5.948501558951e-01;

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

	LinkFileCompatibilitySU3Template(): dim(N,N,N,N), config( dim ), config2( dim )
	{
	}
};

typedef LinkFileCompatibilitySU3Template<4> LinkFileCompatibilitySU3;

#ifdef CULGT_HAVE_LINKFILE_ILDG
TEST_F( LinkFileCompatibilitySU3, CheckILDGPlaquette )
{
	LinkFileILDG<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.ildg" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}
#endif

TEST_F( LinkFileCompatibilitySU3, CheckVogtPlaquette )
{
	LinkFileVogt<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.vogt" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU3, CheckNERSCPlaquette )
{
	LinkFileNERSC<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.nersc" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU3, VogtWriteReadFromNERSC )
{
	LinkFileNERSC<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "lat.sample.l4444.nersc" );
	config2.loadFile( linkfileIn );

	LinkFileVogt<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempILDG2Vogt.vogt" );
	config2.saveFile( linkfileOut );

	LinkFileVogt<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempILDG2Vogt.vogt" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU3, NERSCWriteReadFromVogt )
{
	LinkFileVogt<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "lat.sample.l4444.vogt" );
	config2.loadFile( linkfileIn );


	LinkFileNERSC<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempVogt2NERSC.nersc" );
	config2.saveFile( linkfileOut );

	LinkFileNERSC<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempVogt2NERSC.nersc" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU3, MDPWriteReadFromVogt )
{
	LinkFileVogt<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "lat.sample.l4444.vogt" );
	config2.loadFile( linkfileIn );

	LinkFileMDP<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempVogt2MDP.mdp" );
	config2.saveFile( linkfileOut );

	LinkFileMDP<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempVogt2MDP.mdp" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}


#ifdef CULGT_HAVE_LINKFILE_ILDG
TEST_F( LinkFileCompatibilitySU3, ILDGReadWriteReadWithSameLinkFileWorks )
{
	LinkFileILDG<PATTERNTYPE> linkfileInOut( dim );
	linkfileInOut.setFilename( "lat.sample.l4444.ildg" );
	config2.loadFile( linkfileInOut );

	linkfileInOut.setFilename( "tempILDG2ILDG.ildg" );
	config2.saveFile( linkfileInOut );

	LinkFileILDG<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempILDG2ILDG.ildg" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU3, ILDGReadWriteReadWithDifferentLinkFileThrowsException )
{
	LinkFileILDG<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "lat.sample.l4444.ildg" );
	config2.loadFile( linkfileIn );

	LinkFileILDG<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempILDG2ILDG2.ildg" );

	ASSERT_THROW( config2.saveFile( linkfileOut ), LinkFileException );
}
#endif

typedef LinkFileCompatibilitySU3Template<5> LinkFileCompatibilitySU3WrongSize;

#ifdef CULGT_HAVE_LINKFILE_ILDG
TEST_F( LinkFileCompatibilitySU3WrongSize, ILDGThrowsException )
{
	LinkFileILDG<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.ildg" );

	ASSERT_THROW( config.loadFile( linkfile ), LinkFileException );
}
#endif

TEST_F( LinkFileCompatibilitySU3WrongSize, VogtThrowsException )
{
	LinkFileVogt<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.vogt" );

	ASSERT_THROW( config.loadFile( linkfile ), LinkFileException );
}

TEST_F( LinkFileCompatibilitySU3WrongSize, NERSCThrowsException )
{
	LinkFileNERSC<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "lat.sample.l4444.nersc" );

	ASSERT_THROW( config.loadFile( linkfile ), LinkFileException );
}


template<int Nt, int Nx> class LinkFileCompatibilitySU2Template: public Test
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

	const double plaquetteValue = 0.588338405663378;

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

	LinkFileCompatibilitySU2Template(): dim(Nt,Nx,Nx,Nx), config( dim ), config2( dim )
	{
	}
};

typedef LinkFileCompatibilitySU2Template<16,8> LinkFileCompatibilitySU2;

TEST_F( LinkFileCompatibilitySU2, CheckHirepPlaquette )
{
	LinkFileHirep<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "sample_16x8x8x8_su2.orig.hirep" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU2, CheckVogtPlaquette )
{
	LinkFileVogt<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "sample_16x8x8x8_su2.vogt" );

	config.loadFile( linkfile );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU2, HirepWriteReadFromVogt )
{
	LinkFileVogt<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "sample_16x8x8x8_su2.vogt" );
	config2.loadFile( linkfileIn );

	LinkFileHirep<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempVogt2Hirep.hirep" );
	config2.saveFile( linkfileOut );

	LinkFileHirep<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempVogt2Hirep.hirep" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}

TEST_F( LinkFileCompatibilitySU2, VogtWriteReadFromHirep )
{
	LinkFileHirep<PATTERNTYPE> linkfileIn( dim );
	linkfileIn.setFilename( "sample_16x8x8x8_su2.orig.hirep" );
	config2.loadFile( linkfileIn );


	LinkFileVogt<PATTERNTYPE> linkfileOut( dim );
	linkfileOut.setFilename( "tempHirep2Vogt.vogt" );
	config2.saveFile( linkfileOut );

	LinkFileVogt<PATTERNTYPE> linkfileCheck( dim );
	linkfileCheck.setFilename( "tempHirep2Vogt.vogt" );
	config.loadFile( linkfileCheck );
	config.copyToDevice();

	ASSERT_FLOAT_EQ( plaquetteValue, calcPlaquetteOnConfig() );
}


typedef LinkFileCompatibilitySU2Template<3,8> LinkFileCompatibilitySU2WrongSize;

TEST_F( LinkFileCompatibilitySU2WrongSize, HirepThrowsException )
{
	LinkFileHirep<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "sample_16x8x8x8_su2.orig.hirep" );

	ASSERT_THROW( config.loadFile( linkfile ), LinkFileException );
}

TEST_F( LinkFileCompatibilitySU2WrongSize, VogtThrowsException )
{
	LinkFileVogt<PATTERNTYPE> linkfile( dim );
	linkfile.setFilename( "sample_16x8x8x8_su2.vogt" );

	ASSERT_THROW( config.loadFile( linkfile ), LinkFileException );
}
