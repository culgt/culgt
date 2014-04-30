/**
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef LANDAUGAUGEFIXINGAPP_H_
#define LANDAUGAUGEFIXINGAPP_H_

#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"
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

namespace culgt
{

typedef double REAL;
typedef SU2Vector4<REAL> PARAMTYPE;
typedef SiteIndex<4,FULL_SPLIT> SITE;
typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;
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
		boost::program_options::options_description gaugeOptions("Gaugefixing options");
		gaugeOptions.add_options()
				("fappendix", boost::program_options::value<string>(&fileAppendix)->default_value("gaugefixed_"), "file appendix (append after basename when writing)");

		programOptions->addOption( gaugeOptions );
		landau = new LandauGaugefixing<GLOBALLINK,LOCALLINK>( configuration.getDevicePointer(), dimension );

		landau->orstepsAutoTune<RNG>( programOptions->getSeed() );
	}
private:
	void iterate()
	{
		loadToDevice();
		configuration.copyToDevice();

		GaugeStats stats = landau->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

//		RandomGaugeTrafo<GLOBALLINK,LOCALLINK> random( configuration.getDevicePointer(), dimension );
//		random.randomTrafo( programOptions->getSeed() );

		RunInfo info = landau->orstepsTuned( 1.7, 100 );
		std::cout << "Overrelaxtion: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;

		stats = landau->getGaugeStats();
		std::cout << 1 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		std::cout << "Plaquette: " << plaquette.getPlaquette() << std::endl;

		saveFromDevice( fileAppendix );
	};

	LandauGaugefixing<GLOBALLINK,LOCALLINK>* landau;
	string fileAppendix;
};

} /* namespace culgt */
#endif
