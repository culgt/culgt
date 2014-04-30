/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef COULOMBGAUGEFIXINGAPP_H_
#define COULOMBGAUGEFIXINGAPP_H_

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
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

namespace culgt
{

typedef double REAL;
typedef SU2Vector4<REAL> PARAMTYPE;
typedef SiteIndex<4,TIMESLICE_SPLIT> SITE;
typedef GPUPatternTimesliceParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;
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
		boost::program_options::options_description gaugeOptions("Gaugefixing options");
		gaugeOptions.add_options()
				("fappendix", boost::program_options::value<string>(&fileAppendix)->default_value("gaugefixed_"), "file appendix (append after basename when writing)");

		programOptions->addOption( gaugeOptions );
		coulomb = new CoulombGaugefixing<GLOBALLINKTIMESLICE,LOCALLINK>( configuration.getDevicePointer( 0 ), configuration.getDevicePointer( dim.getDimension(0)-1 ), dim.getDimensionTimeslice(), programOptions->getSeed() );

		coulomb->orstepsAutoTune<RNG>();
		coulomb->sastepsAutoTune<RNG>();
//		coulomb->sastepsAutoTune<RNG>(0.01);
	}
private:
	void iterate()
	{
		loadToDevice();

		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		std::cout << "Plaquette: " << plaquette.getPlaquette() << std::endl;

		coulomb->randomTrafo();

		GaugeStats stats = coulomb->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		RunInfo info;
		info = coulomb->runSimulatedAnnealing( 3, 2, 10 );
		std::cout << "Simulated Annealing: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;

		stats = coulomb->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		info = coulomb->runSimulatedAnnealing( 1.4, 0.01, 1000 );
		std::cout << "Simulated Annealing: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;

		stats = coulomb->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		info = coulomb->runOverrelaxation( 1.7, 1000 );
		std::cout << "Overrelaxtion: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;

		stats = coulomb->getGaugeStats();
		std::cout << 1 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		std::cout << "Plaquette: " << plaquette.getPlaquette() << std::endl;


		GaugeSettings settings;
		settings.setGaugeCopies( 1 );
		settings.setRandomTrafo( true );
		settings.setOrMaxIter( 1000 );
		settings.setOrParameter( 1.8 );

		settings.setSaSteps( 1000 );
		settings.setSaMax( 1.4 );
		settings.setSaMin( 0.1 );

		coulomb->fix( settings );

		stats = coulomb->getGaugeStats();
		std::cout << 1 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;


		saveFromDevice( fileAppendix );
	};

	CoulombGaugefixing<GLOBALLINKTIMESLICE,LOCALLINK>* coulomb;
	string fileAppendix;
};

} /* namespace culgt */
#endif
