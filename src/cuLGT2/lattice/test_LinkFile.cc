


#include "gmock/gmock.h"
#include "LinkFile.h"
#include "testhelper_pattern_stub.h"

using namespace culgt;
using namespace ::testing;

template<typename MemoryConfigurationPattern> class LinkFileStub: public LinkFile<MemoryConfigurationPattern>
{
public:
	MOCK_METHOD0( loadImplementation, void() );
};

TEST( ALinkFile, OpenFileThrowsExceptionIfNoFilenameSet )
{
	LinkFileStub<PatternStub<float> > linkfile;

	ASSERT_THROW( linkfile.openFile(), IOException );
}

TEST( ALinkFile, OpenFileThrowsExceptionIfFileDoesNotExist )
{
	LinkFileStub<PatternStub<float> > linkfile;

	linkfile.setFilename( "fileDoesNotExist.txt" );

	ASSERT_THROW( linkfile.openFile(), FileNotFoundException );
}


TEST( ALinkFile, OpenFileDoesNotThrowExceptionIfFileExists )
{
	LinkFileStub<PatternStub<float> > linkfile;

	linkfile.setFilename( "/dev/random" ); // File that exists

	ASSERT_NO_THROW( linkfile.openFile() );
}

TEST( ALinkFile, OverriddenLoadIsCalledByLoad )
{
	LinkFileStub<PatternStub<float> > linkfile;
	float* U;
	U = new float[10];

	EXPECT_CALL( linkfile, loadImplementation() );

	linkfile.load( U );
}



