/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector2_Real18.h"
#include "application/GaugeConfigurationIteratingApplication.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "lattice/LatticeDimension.h"
#include "observables/PlaquetteAverage.h"
#include "lattice/GlobalLink.h"
#include "gaugefixing/CoulombGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "gaugefixing/RandomGaugeTrafo.h"
#include "gaugefixing/GaugeSettings.h"
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


typedef SiteIndex<4,TIMESLICE_SPLIT> SITE;
typedef GPUPatternTimesliceParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef GlobalLink<PATTERNTYPE,true> GLOBALLINK;
typedef PhiloxWrapper<REAL> RNG;

/*
 *
 */
class CoulombGaugeFixingApp: public GaugeConfigurationIteratingApplication<PATTERNTYPE>
{
public:
	CoulombGaugeFixingApp( const LatticeDimension<PATTERNTYPE::SITETYPE::NDIM> dim, FileIterator fileiterator, ProgramOptions* programOptions )
		: GaugeConfigurationIteratingApplication<PATTERNTYPE>(  dim, fileiterator, programOptions ),
		  plaquette( configuration.getDevicePointer(), dimension )
	{
		boost::program_options::options_description gaugeOptions = settings.getGaugeOptions();
		gaugeOptions.add_options()
				("timeslice", boost::program_options::value<int>(&fixSlice)->default_value(-1), "fix only specific timeslice (-1 = fix all)");
		programOptions->addOption( gaugeOptions );

		coulomb = new CoulombGaugefixing<PATTERNTYPE::TIMESLICE_PATTERNTYPE,LOCALLINK>( configuration.getDevicePointer( 0 ), configuration.getDevicePointer( dim.getDimension(0)-1 ), dim.getDimensionTimeslice(), programOptions->getSeed() );
	}
private:
	GaugeSettings settings;

	void setup()
	{
		coulomb->tune( settings.getTuneFactor() );
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
			if( fixSlice == -1)
				fix();
			else
				fix( fixSlice );
		}
		else
		{
			if( loadToDevice() )
			{
				if( fixSlice == -1)
					fix();
				else
					fix( fixSlice );
				saveFromDevice( settings.getFileAppendix() );
			}
		}
	};

	void fix( int t )
	{
		int tDown = (t == 0)?(dimension.getDimension(0)-1):t-1;
		std::cout << "Timeslice t = " << t << " (" << tDown << ")"<< std::endl;
		coulomb->setTimeslice( configuration.getDevicePointer( t ), configuration.getDevicePointer(tDown) );
		coulomb->fix( settings );
	}

	void fix()
	{
		std::cout << "Plaquette before: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;

		RunInfo info;

		for( int t = 0; t < dimension.getDimension(0); t++ )
		{
			fix(t);
		}

		std::cout << "Plaquette after: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;
	}

	CoulombGaugefixing<PATTERNTYPE::TIMESLICE_PATTERNTYPE,LOCALLINK>* coulomb;
	PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette;
	int fixSlice;
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
	CoulombGaugeFixingApp::main<CoulombGaugeFixingApp>( argc, argv );
}
