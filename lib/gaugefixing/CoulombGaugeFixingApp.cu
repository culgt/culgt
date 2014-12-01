/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "application/GaugeConfigurationIteratingApplication.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "lattice/LinkFile.h"
#include "lattice/LatticeDimension.h"
#include "lattice/LinkFileVogt.h"
#include "observables/PlaquetteAverage.h"
#include "lattice/GlobalLink.h"
#include "gaugefixing/CoulombGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "gaugefixing/RandomGaugeTrafo.h"
#include "gaugefixing/GaugeSettings.h"
//#include "lattice/parameterization_types/SUNRealFull.h"

namespace culgt
{
#ifdef DOUBLEPRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

//typedef SU3Vector4<REAL> PARAMTYPE;
//typedef LocalLink<SUNRealFull<3,REAL> > LOCALLINK;

typedef SU2Vector4<REAL> PARAMTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;


typedef SiteIndex<4,TIMESLICE_SPLIT> SITE;
typedef GPUPatternTimesliceParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef GlobalLink<PATTERNTYPE,true> GLOBALLINK;
typedef GlobalLink<PATTERNTYPE::TIMESLICE_PATTERNTYPE,true> GLOBALLINKTIMESLICE;
typedef PhiloxWrapper<REAL> RNG;

/*
 *
 */
class CoulombGaugeFixingApp: public GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >
{
public:
	CoulombGaugeFixingApp( const LatticeDimension<PATTERNTYPE::SITETYPE::Ndim> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >(  dim, fileiterator, programOptions )
	{
		programOptions->addOption( settings.getGaugeOptions() );

		boost::program_options::options_description gaugeOptions;
		gaugeOptions.add_options()
				("fappendix", boost::program_options::value<string>(&fileAppendix)->default_value("gaugefixed_"), "file appendix (append after basename when writing)");

		programOptions->addOption( gaugeOptions );


		coulomb = new CoulombGaugefixing<GLOBALLINKTIMESLICE,LOCALLINK>( configuration.getDevicePointer( 0 ), configuration.getDevicePointer( dim.getDimension(0)-1 ), dim.getDimensionTimeslice(), programOptions->getSeed() );
	}
private:
	GaugeSettings settings;

	void setup()
	{
		coulomb->orstepsAutoTune<RNG>(1.5, 200);
		coulomb->cornellAutoTune<RNG>(.5, 200);
		coulomb->sastepsAutoTune<RNG>(.5, 200);
		coulomb->microcanonicalAutoTune<RNG>( 200 );
	}

	void teardown()
	{
	}

	void iterate()
	{
		loadToDevice();

//		configuration.setColdOnDevice();

		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		std::cout << "Plaquette: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;

		int t = 0;
		int tDown = (t == 0)?(dimension.getDimension(0)-1):t-1;
		std::cout << "Timeslice t = " << t << " (" << tDown << ")"<< std::endl;
		coulomb->setTimeslice( configuration.getDevicePointer( t ), configuration.getDevicePointer(tDown) );

		GaugeStats stats = coulomb->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;
//		stats = coulomb->getGaugeStats( GAUGEFIELD_LOGARITHMIC );
//		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		coulomb->randomTrafo();

		coulomb->fix( settings );
//		for( int j = 0; j < 10000; j++ )
//		{
//			int iter = 100;
//			for( int i = 0; i < iter; i++ )
//			{
//				coulomb->runCornell( .1, 5 );
////				coulomb->runOverrelaxation( 1.5 );
//				CUDA_LAST_ERROR( "Cornell ");
//			}
//				coulomb->reproject();
//
//
//			stats = coulomb->getGaugeStats();
//			std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;
//			if( stats.getPrecision() < settings.getPrecision() ) break;
//			stats = coulomb->getGaugeStats( GAUGEFIELD_LOGARITHMIC );
//			std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;
//		}
		std::cout << "Plaquette: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;
//		exit( 1 );

		RunInfo info;

//		for( int t = 0; t < dimension.getDimension(0); t++ )
//		{
//			int tDown = (t == 0)?(dimension.getDimension(0)-1):t-1;
//			std::cout << "Timeslice t = " << t << " (" << tDown << ")"<< std::endl;
//			coulomb->setTimeslice( configuration.getDevicePointer( t ), configuration.getDevicePointer(tDown) );
////			coulomb->randomTrafo();
//			coulomb->fix( settings );
//		}

		std::cout << "Plaquette: " << std::setprecision(12) << plaquette.getPlaquette() << std::endl;

		saveFromDevice( fileAppendix );
	};

	CoulombGaugefixing<GLOBALLINKTIMESLICE,LOCALLINK>* coulomb;
	string fileAppendix;
};

} /* namespace culgt */


using namespace culgt;

int main( const int argc, const char* argv[] )
{
	CoulombGaugeFixingApp::main<CoulombGaugeFixingApp>( argc, argv );
}
