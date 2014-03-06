/**
 * test_LinkFileWithPattern.cc
 *
 *  Created on: Mar 6, 2014
 *      Author: vogt
 */

#include "gmock/gmock.h"
#include "LinkFileVogt.h"
#include "configuration_patterns/GPUPattern.h"
#include "su3/SU3Real12.h"
#include "../cuLGT1legacy/SiteIndex.hxx"

using namespace culgt;
using namespace ::testing;


/**
 * Checks that the first entry (real part of 1/1 entry at site 0 in mu=0 direction)
 */
//TEST( ALinkFileVogtWithPatternLoadingSampleConfiguration, FirstDataEntryAfterHeaderIsCorrect )
//{
//	float* U = new float[8*4*4*4*4*8];
//	LinkFileVogt<GPUPattern<SiteIndex<4,NO_SPLIT>, SU3Real12<float> >, float > linkfile;
//
//	linkfile.setPointerToU(U);
//	linkfile.loadHeader();
//
//	linkfile.loadBody();
//
//	ASSERT_EQ( -1.966492e-01, U[0] );
//}

