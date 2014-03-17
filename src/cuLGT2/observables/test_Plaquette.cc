/**
 * test_Plaquette.cc
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "../lattice/GaugeConfigurationHelper.h"
#include "../lattice/configuration_patterns/StandardPattern.h"
#include "../lattice/parameterization_types/SUNRealFull.h"
#include "../cuLGT1legacy/SiteCoord.hxx"
#include "Plaquette.h"


using namespace culgt;
using namespace ::testing;


class APlaquette: public Test
{
public:
	static const int Ndim = 2;
	static const int N = 2;
	static const int mu = 0;
	static const int nu = 1;
	typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>, SUNRealFull<3,float> > MyPattern;

	float U[N*N*Ndim*SUNRealFull<3,float>::SIZE];
	const LatticeDimension<Ndim> dim;
	Plaquette<MyPattern> plaquette;
	MyPattern::SITETYPE site;

	APlaquette() : dim(N,N), plaquette( U, dim ), site( dim )
	{
		GaugeConfigurationHelper<MyPattern>::setCold( U, dim );
		site.setLatticeIndex(0);
	}

};

TEST_F( APlaquette, GetStapleForCold2x2LatticeReturnsIdentity )
{
	ASSERT_TRUE( LocalLink<MyPattern::PARAMTYPE>::getIdentity() == plaquette.getStaple( site, mu, nu ) );
}

TEST_F( APlaquette, GetPlaquetteForCold2x2LatticeReturnsIdentity )
{
	ASSERT_TRUE( LocalLink<MyPattern::PARAMTYPE>::getIdentity() == plaquette.getPlaquette( site, mu, nu ) );
}

TEST_F( APlaquette, GetReTracePlaquetteNormalizedForCold2x2LatticeReturns1 )
{
	ASSERT_FLOAT_EQ( 1.0, plaquette.getReTracePlaquetteNormalized( site, mu, nu ) );
}


