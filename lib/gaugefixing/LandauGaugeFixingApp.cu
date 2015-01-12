/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"
#include "application/GaugeConfigurationIteratingApplication.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "lattice/LinkFile.h"
#include "lattice/LatticeDimension.h"
#include "lattice/LinkFileVogt.h"
#include "observables/PlaquetteAverage.h"
#include "lattice/GlobalLink.h"
#include "gaugefixing/LandauGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "lattice/parameterization_types/SU3Vector4.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/SU3Real12.h"

namespace culgt
{

#ifdef DOUBLEPRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

#ifdef CULGT_SU2
typedef SU2Vector4<REAL> PARAMTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;
#else
typedef SU3Vector4<REAL> PARAMTYPE;
//typedef SU3Real12<REAL> PARAMTYPE;
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

/*
 *
 */
class LandauGaugeFixingApp: public GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >
{
public:
	LandauGaugeFixingApp( const LatticeDimension<PATTERNTYPE::SITETYPE::Ndim> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >(  dim, fileiterator, programOptions )
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
	}

	void teardown()
	{
	}

	void iterate()
	{
		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		if( sethot )
		{
			configuration.setHotOnDevice<RNG>( programOptions->getSeed(), RNG::getNextCounter());
			CUDA_LAST_ERROR( "setHotOnDevice ");
		}
		else
		{
			loadToDevice();
		}

		std::cout << "Plaquette before: " << plaquette.getPlaquette() << std::endl;
		landau->fix( settings );
		std::cout << "Plaquette after: " << plaquette.getPlaquette() << std::endl;

		saveFromDevice( fileAppendix );
	};

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
	LandauGaugeFixingApp::main<LandauGaugeFixingApp>( argc, argv );
}