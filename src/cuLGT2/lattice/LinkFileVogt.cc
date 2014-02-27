/**
 * LinkFile.h
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef LINKFILEVOGT_H_
#define LINKFILEVOGT_H_

#include "LinkFile.h"
#include <string>
#include <fstream>

namespace culgt
{

template<typename MemoryConfigurationPattern, typename TFloatFile> class LinkFileVogt: LinkFile<MemoryConfigurationPattern>
{

public:
	virtual void loadImplementation() override;

private:
	void loadHeader()
	{
	};
};

template<typename MemoryConfigurationPattern, typename TFloatFile > void LinkFileVogt<MemoryConfigurationPattern, TFloatFile>::loadImplementation()
{
}


}

#endif /* LINKFILE_H_ */


#include "gmock/gmock.h"
#include "testhelper_pattern_stub.h"

using namespace culgt;
using namespace ::testing;


//TEST( ALinkFileVogt, OpenFileThrowsExceptionIfFilenameNotSet )
//{
//	LinkFileVogt<PatternStub<float>,float> linkfile;
//	linkfile.load( U );
//	FAIL();
//}


