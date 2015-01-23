/*
 * MILCcuLGT.cpp
 *
 *  Created on: Aug 18, 2014
 *      Author: vogt
 */

#include "MILCcuLGT.h"
#include "lattice/GaugeConfiguration.h"
#include <iostream>
#include "observables/WilsonLoopAverage.h"
#include "lattice/parameterization_types/SU3Vector4.h"
#include "cuLGT1legacy/SiteIndex.hxx"
#include "lattice/configuration_patterns/GPUPatternParityPriority.h"
#include "lattice/GlobalLink.h"
#include "lattice/LocalLink.h"
#include "interfaces/milc/MILCConverter.h"
#include  "lattice/parameterization_types/ParameterizationMediatorSU3_Vector4_Real18.h"
#include "gaugefixing/LandauGaugeFixing.h"
#include "gaugefixing/GaugeSettings.h"
#include "util/rng/PhiloxWrapper.h"

#include "su3.h" // milc
#include "lattice.h"

using culgt::GaugeConfiguration;
using culgt::SU3Vector4;
using culgt::SiteIndex;
using culgt::GPUPatternParityPriority;
using culgt::GlobalLink;
using culgt::LocalLink;
using culgt::MILCConverter;
using culgt::LandauGaugefixing;
using culgt::GaugeSettings;
using culgt::PhiloxWrapper;

typedef SU3Vector4<REAL> PARAMTYPE;
typedef SiteIndex<4,FULL_SPLIT> SITE;
typedef GPUPatternParityPriority<SITE,PARAMTYPE> PATTERNTYPE;
typedef GlobalLink<PATTERNTYPE,true> GLOBALLINK;
typedef LocalLink<SUNRealFull<3,REAL> > LOCALLINK;
typedef PhiloxWrapper<REAL> RNG;

LandauGaugefixing<PATTERNTYPE,LOCALLINK>* landau;
GaugeConfiguration<PATTERNTYPE>* config;

void cuLGTinitLandau( int nx, int ny, int nz, int nt  )
{
	LatticeDimension<4> dim( nt, nx, ny, nz );
	config = new GaugeConfiguration<PATTERNTYPE>(dim);
	config->allocateMemory();
	landau = new LandauGaugefixing<PATTERNTYPE,LOCALLINK>( config->getDevicePointer(), dim, 1235 );
	landau->orstepsAutoTune<RNG>(1.5, 200);
	landau->microcanonicalAutoTune( 200 );
	landau->sastepsAutoTune<RNG>(1., 20);
}

void cuLGTfixLandau( int nx, int ny, int nz, int nt )
{
	MILCConverter<PATTERNTYPE, REAL>::convertFromMILC( config->getHostPointer(), nx, ny, nz, nt );

	config->copyToDevice();

	GaugeSettings settings;
	settings.setOrMaxIter( 6000 );
	settings.setSaSteps( 200 );
	settings.setSaMin( 0.01 );
	settings.setSaMax( 1.4 );
	settings.setMicroiter( 3 );
	settings.setReproject( 100 );
	settings.setOrParameter( 1.5 );
	settings.setGaugeCopies( 1 );
	settings.setPrecision( 1e-6 );
	settings.setPrintStats( true );
	settings.setCheckPrecision( 100 );

	landau->fix( settings );

	config->copyToHost();

	MILCConverter<PATTERNTYPE, REAL>::convertToMILC( config->getHostPointer(), nx, ny, nz, nt );
}
