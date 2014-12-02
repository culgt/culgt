


#include "gmock/gmock.h"
#include "lattice/LinkFile.h"
#include "testhelper_pattern_stub.h"

using namespace culgt;
using namespace ::testing;

template<typename MemoryConfigurationPattern> class LinkFileStub: public LinkFile<MemoryConfigurationPattern>
{
public:
	MOCK_METHOD0( loadImplementation, void() );
	MOCK_METHOD0( saveImplementation, void() );
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



class ALinkFileWithDestination: public Test
{
public:
	LinkFileStub<PatternStub<float> > linkfile;
	float* U;

	void SetUp()
	{
		U = new float[1];
		linkfile.setFilename( "/dev/random" ); // File that exists
	}
};

TEST_F( ALinkFileWithDestination, LoadThrowsIOExceptionBecauseItCallsOpen )
{
	LinkFileStub<PatternStub<float> > linkfile;

	ASSERT_THROW( linkfile.load( U ), IOException );
}

TEST_F( ALinkFileWithDestination, OverriddenLoadIsCalledByLoad )
{
	EXPECT_CALL( linkfile, loadImplementation() );
	linkfile.load( U );
}

TEST_F( ALinkFileWithDestination, OverriddenSaveIsCalledByLoad )
{
	EXPECT_CALL( linkfile, saveImplementation() );
	linkfile.save( U );
}



