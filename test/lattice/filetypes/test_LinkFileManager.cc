/**
 * test_LinkFileWithPattern.cc
 *
 *  Created on: Mar 6, 2014
 *      Author: vogt
 */

#include "gmock/gmock.h"
#include "lattice/filetypes/LinkFileManager.h"
#include "../testhelper_pattern_stub.h"

using namespace culgt;
using namespace ::testing;


TEST( ALinkFileManager, AllocationForDifferentFileTypes )
{
	typedef PatternStub<float,4,2> MEMORYPATTERN;
	LatticeDimension<4> dim;

	LinkFileManager<MEMORYPATTERN> filemanagerVogt( VOGT, dim, STANDARD );
	LinkFileManager<MEMORYPATTERN> filemanagerHeaderonly( HEADERONLY, dim, STANDARD );
	LinkFileManager<MEMORYPATTERN> filemanagerHirep( HIREP, dim, STANDARD );
	LinkFileManager<MEMORYPATTERN> filemanagerILDG( ILDG, dim, STANDARD );
}
