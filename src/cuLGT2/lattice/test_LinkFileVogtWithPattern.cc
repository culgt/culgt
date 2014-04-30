/**
 * test_LinkFileWithPattern.cc
 *
 *  Created on: Mar 6, 2014
 *      Author: vogt
 */

#include "gmock/gmock.h"
#include "LinkFileVogt.h"
#include "configuration_patterns/GPUPattern.h"
#include "parameterization_types/ParameterizationMediatorSU2_Real4_Real8.h"
#include "../cuLGT1legacy/SiteIndex.hxx"
#include "LatticeDimension.h"
#include "../common/culgt_typedefs.h"

using namespace culgt;
using namespace ::testing;


class ALinkFileVogtWithPatternLoadingSampleBody: public Test
{
public:
	float* U = new float[8*4*4*4*4*8];
	LinkFileVogt<GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> >, float >* linkfile;

	void SetUp()
	{
		LatticeDimension<4> dim(8,4,4,4);
		linkfile = new LinkFileVogt<GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> >, float >( dim );
		linkfile->setFilename( "test_configSU2N4T8SP.vogt" );
		linkfile->openFile();

		linkfile->setPointerToU(U);
		linkfile->loadHeader();
		linkfile->verify();
	}
};

/**
 * Checks that the first entry (real part of 1/1 entry at site 0 in mu=0 direction)
 */
TEST_F( ALinkFileVogtWithPatternLoadingSampleBody, FirstDataEntryAfterHeaderIsCorrect )
{
	linkfile->loadBody();

	ASSERT_FLOAT_EQ( -4.711566e-01, U[0] );
}

//class ALinkFileVogtWithPatternLoadingSample: public Test
//{
//public:
//
//	float* U = new float[8*4*4*4*4*8];
//	LinkFileVogt<GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> >, float >* linkfile;
//
//	void SetUp()
//	{
//		LatticeDimension<4> dim(8,4,4,4);
//		linkfile = new LinkFileVogt<GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> >, float >( dim );
//		linkfile->setFilename( "../test_configSU2N4T8SP.vogt" );
//	}
//};

TEST( ALinkFileVogtWithPatternLoadingSample, LoadViaAbstractLoadReadsCorrectFirstEntry )
{
	typedef GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> > MyPattern;
	float* U = new float[8*4*4*4*4*8];
	LatticeDimension<4> dim(8,4,4,4);
	LinkFile<MyPattern>* aLinkFileVogt;
	aLinkFileVogt = new LinkFileVogt<MyPattern,float>( dim );
	aLinkFileVogt->setFilename( "test_configSU2N4T8SP.vogt" );

	aLinkFileVogt->load( U );

	ASSERT_FLOAT_EQ( -4.711566e-01, U[0] );
}

TEST( ALinkFileVogtWithPatternSavingSample, LoadAfterSaveWorks )
{
	const float someValue = 54.2345;

	typedef GPUPattern<SiteIndex<4,NO_SPLIT>, SU2Real4<float> > MyPattern;
	float* U = new float[8*4*4*4*4*8];
	LatticeDimension<4> dim(8,4,4,4);
	LinkFile<MyPattern>* aLinkFileVogt;
	aLinkFileVogt = new LinkFileVogt<MyPattern,float>( dim );

	U[0] = someValue;

	aLinkFileVogt->setFilename( "test_save.vogt" );

	aLinkFileVogt->save( U );

	U[0] = 0;
	aLinkFileVogt->load( U );


	ASSERT_FLOAT_EQ( someValue, U[0] );
}


