/*
 * testhelper_linkfile_stub.h
 *
 *  Created on: Jul 3, 2015
 *      Author: vogt
 */

#ifndef TESTHELPER_LINKFILE_STUB_H_
#define TESTHELPER_LINKFILE_STUB_H_

#include "gmock/gmock.h"
#include <string>
#include "../testhelper_pattern_stub.h"
#include "lattice/filetypes/LinkFile.h"

using std::string;
using culgt::LinkFile;

template<typename MemoryConfigurationPattern> class LinkFileStub: public LinkFile<MemoryConfigurationPattern>
{
public:
	MOCK_METHOD0( loadImplementation, void() );
	MOCK_METHOD0( saveImplementation, void() );
	MOCK_METHOD0( getPreferredExtension, string() );
};

typedef LinkFileStub<PatternStub<float> > LinkFileStubDefault;

#endif /* TESTHELPER_LINKFILE_STUB_H_ */
