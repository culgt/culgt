/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector2_Real18.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"
#include "application/GaugeConfigurationIteratingApplication.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "lattice/LinkFile.h"
#include "lattice/LatticeDimension.h"
#include "lattice/LinkFileVogt.h"
#include "lattice/LinkFileHirep.h"
#include "lattice/LinkFileHeaderOnly.h"
#include "observables/PlaquetteAverage.h"
#include "lattice/GlobalLink.h"
#include "gaugefixing/LandauGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "version.h"

#if __cplusplus >= 201103L
// this fixes a strange error in boost/lexical_cast.hpp where it needs std::pow( double, int )
#include <cmath>
double std::pow( double d, int i )
{
	return std::pow( d, (double)i );
}
#endif

namespace culgt
{

#ifdef DOUBLEPRECISION
typedef double REAL;
#else
typedef float REAL;
#endif


#if CULGT_SUN == 2

typedef SU2Vector4<REAL> PARAMTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;

#else

#ifdef DOUBLEPRECISION
typedef SU3Vector2<REAL> PARAMTYPE;
#else
typedef SU3Vector4<REAL> PARAMTYPE;
#endif
typedef LocalLink<SUNRealFull<3,REAL> > LOCALLINK;

#endif


#ifdef CULGT_USE_TIMESLICE_PATTERN
typedef SiteIndex<4,TIMESLICE_SPLIT> SITE;
typedef GPUPatternTimesliceParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
#else
typedef SiteIndex<4,FULL_SPLIT> SITE;
typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
#endif

typedef GlobalLink<PATTERNTYPE,true> GLOBALLINK;
typedef PhiloxWrapper<REAL> RNG;

#ifdef CULGT_FILETYPE_VOGT
typedef LinkFileVogt<PATTERNTYPE,REAL> FILETYPE;
#elif CULGT_FILETYPE_HIREP
typedef LinkFileHirep<PATTERNTYPE,REAL> FILETYPE;
#else
typedef LinkFileHeaderOnly<PATTERNTYPE,REAL> FILETYPE;
#endif

#ifdef CULGT_FILETYPE_VOGT_OUT
typedef LinkFileVogt<PATTERNTYPE,REAL> FILETYPEOUT;
#elif CULGT_FILETYPE_HIREP_OUT
typedef LinkFileHirep<PATTERNTYPE,REAL> FILETYPEOUT;
#elif CULGT_FILETYPE_HEADERONLY_OUT
typedef LinkFileHeaderOnly<PATTERNTYPE,REAL> FILETYPEOUT;
#else
typedef FILETYPE FILETYPEOUT;
#endif

/*
 *
 */
class LandauGaugeFixingApp: public GaugeConfigurationIteratingApplication<PATTERNTYPE,FILETYPE,FILETYPEOUT>
{
public:
	LandauGaugeFixingApp( const LatticeDimension<PATTERNTYPE::SITETYPE::Ndim> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : GaugeConfigurationIteratingApplication<PATTERNTYPE,FILETYPE,FILETYPEOUT>(  dim, fileiterator, programOptions ), plaquette( configuration.getDevicePointer(), dimension )
	{
		programOptions->addOption( settings.getGaugeOptions() );

		boost::program_options::options_description gaugeOptions("Gaugefixing options");
		gaugeOptions.add_options()
				("sethot", boost::program_options::value<bool>(&sethot)->default_value(false), "start from a random gauge field")
				("fappendix", boost::program_options::value<string>(&fileAppendix)->default_value("gaugefixed_"), "file appendix (append after basename when writing)");

		programOptions->addOption( gaugeOptions );
		landau = new LandauGaugefixing<PATTERNTYPE,LOCALLINK>( configuration.getDevicePointer(), dimension, programOptions->getSeed() );
	}
private:
	GaugeSettings settings;

	void setup()
	{
		landau->orstepsAutoTune<RNG>(1.5, 50);
		landau->microcanonicalAutoTune<RNG>( 50 );
		landau->sastepsAutoTune<RNG>(1., 20);
	}

	void teardown()
	{
	}

	void iterate()
	{
		if( sethot )
		{
			std::cout << "Using a hot (random) lattice" << std::endl;
			configuration.setHotOnDevice<RNG>( programOptions->getSeed(), RNG::getNextCounter());
			CUDA_LAST_ERROR( "setHotOnDevice ");
			fix();
		}
		else
		{
			if( loadToDevice() )
			{
				fix();
				saveFromDevice( fileAppendix );
			}
		}
	};

	void fix()
	{
		std::cout << "Plaquette before: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;
		landau->fix( settings );
		std::cout << "Plaquette after: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;
	}

	PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette;
	LandauGaugefixing<PATTERNTYPE,LOCALLINK>* landau;
	string fileAppendix;
	bool sethot;
};

} /* namespace culgt */


using namespace culgt;

#if BOOST_VERSION < 105300
int main( int argc, char* argv[] )
#else
int main( const int argc, const char* argv[] )
#endif
{
	std::cout << "cuLGT Version " << CULGT_VERSION << std::endl;
	LandauGaugeFixingApp::main<LandauGaugeFixingApp>( argc, argv );
}
