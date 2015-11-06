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
#include "lattice/site_indexing/SiteIndex.h"
#include "lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "lattice/LatticeDimension.h"
#include "observables/PlaquetteAverage.h"
#include "lattice/GlobalLink.h"
#include "gaugefixing/LandauGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "version.h"
#include "common/compiler.h"

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


/*
 *
 */
class LandauGaugeFixingApp: public GaugeConfigurationIteratingApplication<PATTERNTYPE>
{
public:
	LandauGaugeFixingApp( const LatticeDimension<PATTERNTYPE::SITETYPE::NDIM> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : GaugeConfigurationIteratingApplication<PATTERNTYPE>(  dim, fileiterator, programOptions ), plaquette( configuration.getDevicePointer(), dimension )
	{
		programOptions->addOption( settings.getGaugeOptions() );
		landau = new LandauGaugefixing<PATTERNTYPE,LOCALLINK>( configuration.getDevicePointer(), dimension, programOptions->getSeed() );
	}
private:
	GaugeSettings settings;
	ofstream perfout;

	void setup()
	{
		landau->tune( settings.getTuneFactor() );
	}

	void teardown()
	{
	}

	void iterate()
	{
		if( settings.isSethot() )
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
				saveFromDevice( settings.getFileAppendix() );
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
	std::cout << "Compiled with " << CULGT_NVCC_COMPILER << " (" << CULGT_HOST_COMPILER << ") using boost version " << BOOST_VERSION << std::endl;
	LandauGaugeFixingApp::main<LandauGaugeFixingApp>( argc, argv );
}
