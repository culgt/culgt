/**
 * CoulombGaugeFixingApp2.h
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef COULOMBGAUGEFIXINGAPP2_H_
#define COULOMBGAUGEFIXINGAPP2_H_

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

namespace culgt
{

typedef float REAL;
typedef SU2Vector4<REAL> PARAMTYPE;
typedef SiteIndex<4,FULL_SPLIT> SITE;
typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef LocalLink<SU2Vector4<REAL> > LOCALLINK;
typedef GlobalLink<PATTERNTYPE,false> GLOBALLINK;

/*
 *
 */
class CoulombGaugeFixingApp2: public GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >
{
public:
	CoulombGaugeFixingApp2( const LatticeDimension<PATTERNTYPE::SITETYPE::Ndim> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : GaugeConfigurationIteratingApplication<PATTERNTYPE,LinkFileVogt<PATTERNTYPE,REAL> >(  dim, fileiterator, programOptions )
	{
		landau = new LandauGaugefixing<GLOBALLINK,LOCALLINK>( configuration.getDevicePointer(), dimension );
	}
private:
	void iterate()
	{
		load();
		configuration.copyToDevice();

		GaugeStats stats = landau->getGaugeStats();
		std::cout << 0 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		landau->orsteps<EIGHT_THREAD_PER_SITE,64,5>( 1.7, 1 );
		RunInfo info = landau->orsteps<EIGHT_THREAD_PER_SITE,64,5>( 1.7, 100 );
		std::cout << "Overrelaxtion: " << info.getGflops() << " GFlops at " <<  info.getThroughput() << " GByte/s memory throughput." << std::endl;

		stats = landau->getGaugeStats();
		std::cout << 1 << " \t" << stats.getGff() << " \t" << stats.getPrecision() << std::endl;

		PlaquetteAverage<PATTERNTYPE,LOCALLINK> plaquette( configuration.getDevicePointer(), dimension );
		std::cout << "Plaquette: " << plaquette.getPlaquette() << std::endl;
	};

	LandauGaugefixing<GLOBALLINK,LOCALLINK>* landau;
};

} /* namespace culgt */
#endif /* COULOMBGAUGEFIXINGAPP2_H_ */
