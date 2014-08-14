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
#include "gaugefixing/LandauGaugeFixing.h"
#include "util/rng/PhiloxWrapper.h"
#include "lattice/parameterization_types/SU3Vector4.h"
#include "lattice/parameterization_types/SUNRealFull.h"

namespace culgt
{

#ifdef DOUBLEPRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

//typedef SU2Vector4<REAL> PARAMTYPE;
//typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;

typedef SU3Vector4<REAL> PARAMTYPE;
typedef LocalLink<SUNRealFull<3,REAL> > LOCALLINK;


typedef SiteIndex<4,FULL_SPLIT> SITE;
typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
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
		landau = new LandauGaugefixing<GLOBALLINK,LOCALLINK>( configuration.getDevicePointer(), dimension, programOptions->getSeed() );
	}
private:
	GaugeSettings settings;

	void setup()
	{
		landau->orstepsAutoTune<RNG>(1.5, 200);
	}

	void teardown()
	{
	}

	void iterate()
	{
		if( sethot )
		{
			configuration.setHotOnDevice<RNG>( programOptions->getSeed(), RNG::getNextCounter());
			CUDA_LAST_ERROR( "setHotOnDevice ");
		}
		else
		{
			loadToDevice();
		}
//		configuration.copyToDevice();

		landau->fix( settings );

//		GaugeStats stats = landau->getGaugeStats();
//		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;
//
////		RandomGaugeTrafo<GLOBALLINK,LOCALLINK> random( configuration.getDevicePointer(), dimension );
////		random.randomTrafo( programOptions->getSeed() );
//
//		RunInfo info = landau->orstepsTuned( 1.7, 100 );
//		std::cout << "Overrelaxtion: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;
//
//		stats = landau->getGaugeStats();
//		std::cout << 1 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		std::cout << "Plaquette: " << plaquette.getPlaquette() << std::endl;

		saveFromDevice( fileAppendix );
	};

	LandauGaugefixing<GLOBALLINK,LOCALLINK>* landau;
	string fileAppendix;
	bool sethot;
};

} /* namespace culgt */


using namespace culgt;

int main( const int argc, const char* argv[] )
{
	LandauGaugeFixingApp::main<LandauGaugeFixingApp>( argc, argv );
}
