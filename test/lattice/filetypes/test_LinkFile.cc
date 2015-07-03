


#include "gmock/gmock.h"
#include "lattice/filetypes/LinkFile.h"
#include "../testhelper_pattern_stub.h"
#include "testhelper_linkfile_stub.h"

using namespace culgt;
using namespace ::testing;


TEST( ALinkFile, OpenFileThrowsExceptionIfNoFilenameSet )
{
	LinkFileStubDefault linkfile;

	ASSERT_THROW( linkfile.openFile(), IOException );
}

TEST( ALinkFile, OpenFileThrowsExceptionIfFileDoesNotExist )
{
	LinkFileStubDefault linkfile;

	linkfile.setFilename( "fileDoesNotExist.txt" );

	ASSERT_THROW( linkfile.openFile(), FileNotFoundException );
}


TEST( ALinkFile, OpenFileDoesNotThrowExceptionIfFileExists )
{
	LinkFileStubDefault linkfile;

	linkfile.setFilename( "/dev/random" ); // File that exists

	ASSERT_NO_THROW( linkfile.openFile() );
}



class ALinkFileWithDestination: public Test
{
public:
	LinkFileStubDefault linkfile;
	float* U;

	void SetUp()
	{
		U = new float[1];
		linkfile.setFilename( "/dev/random" ); // File that exists
	}
};

TEST_F( ALinkFileWithDestination, LoadThrowsIOExceptionBecauseItCallsOpen )
{
	LinkFileStubDefault linkfile;

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



