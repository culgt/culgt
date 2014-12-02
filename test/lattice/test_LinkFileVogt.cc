
#include "gmock/gmock.h"
#include "testhelper_pattern_stub.h"
#include "lattice/LinkFileVogt.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"

using namespace culgt;
using namespace ::testing;

/**
 * The sample configuration is 8x4^3 lattice SU(2) lattice in Single Precision with first entry -4.711566e-01
 */
class ALinkFileVogtWithSampleConfiguration: public Test
{
public:
	LinkFileVogt<PatternStub<float,4,2>,float> linkfile;

	void SetUp()
	{
		linkfile.setFilename("test_configSU2N4T8SP.vogt");
		linkfile.openFile();
	}
	void TearDown()
	{
		linkfile.closeFile();
	}
};

TEST_F( ALinkFileVogtWithSampleConfiguration, LoadHeaderReadsNdim )
{
	linkfile.loadHeader();

	ASSERT_EQ( 4, linkfile.getNdim() );
}

TEST_F( ALinkFileVogtWithSampleConfiguration, LoadHeaderReadsNc )
{
	linkfile.loadHeader();

	ASSERT_EQ( 2, linkfile.getNc() );
}

TEST_F( ALinkFileVogtWithSampleConfiguration, LoadHeaderReadsLatticeSize )
{
	linkfile.loadHeader();

	ASSERT_EQ( 8, linkfile.getSize(0) );
	ASSERT_EQ( 4, linkfile.getSize(3) );
}

TEST_F( ALinkFileVogtWithSampleConfiguration, LoadHeaderReadsSizeOfReal)
{
	linkfile.loadHeader();

	ASSERT_EQ( sizeof( float ), (size_t)linkfile.getSizeOfReal() );
}

TEST_F( ALinkFileVogtWithSampleConfiguration, LoadFirstLinkHasCorrectFirstEntry)
{
	const float firstEntry = -4.711566e-01;
	linkfile.loadHeader();
	LocalLink<SUNRealFull<2,float> >link;
	link = linkfile.getNextLink();

	ASSERT_EQ( firstEntry, link.get(0) );
}

template<typename LinkFileType> void readSampleFileHeader( LinkFileType& linkfile )
{
	linkfile.setFilename("test_configSU2N4T8SP.vogt");
	linkfile.openFile();
	linkfile.loadHeader();
}

TEST( ALinkFileVogtWithWrongSettings, VerifyThrowsExceptionIfWrongNdim )
{
	const int Nc = 2;
	const int WrongNDIM = 3;
	LinkFileVogt<PatternStub<float,WrongNDIM,Nc>,float> linkfile;

	readSampleFileHeader( linkfile );

	try
	{
		linkfile.verify();
		FAIL();
	}
	catch( LinkFileVogtException e )
	{
		ASSERT_THAT( e.what(), StartsWith("Wrong lattice dimension") );
	}
}

TEST( ALinkFileVogtWithWrongSettings, VerifyThrowsExceptionIfWrongNc )
{
	const int NDIM = 4;
	const int WrongNc = 3;
	LinkFileVogt<PatternStub<float,NDIM, WrongNc>,float> linkfile;

	readSampleFileHeader( linkfile );

	try
	{
		linkfile.verify();
		FAIL();
	}
	catch( LinkFileVogtException e )
	{
		ASSERT_THAT( e.what(), StartsWith("Wrong gauge group") );
	}
}

TEST( ALinkFileVogtWithWrongSettings, VerifyThrowsExceptionIfWrongSizeOfReal)
{
	const int NDIM = 4;
	const int Nc = 2;
	LinkFileVogt<PatternStub<double,NDIM,Nc>,double> linkfile;

	readSampleFileHeader( linkfile );

	try
	{
		linkfile.verify();
		FAIL();
	}
	catch( LinkFileVogtException e )
	{
		ASSERT_THAT( e.what(), StartsWith("Wrong size of real") );
	}
}

TEST( ALinkFileVogtWithWrongSettings, VerifyThrowsExceptionIfWrongLatticeSize)
{
	const int NDIM = 4;
	const int Nc = 2;
	const int size[4] = {7,4,4,4};
	LinkFileVogt<PatternStub<float,NDIM,Nc>,float> linkfile( size );

	readSampleFileHeader( linkfile );

	try
	{
		linkfile.verify();
		FAIL();
	}
	catch( LinkFileVogtException e )
	{
		ASSERT_THAT( e.what(), StartsWith("Wrong lattice size in 0 direction") );
	}
}

TEST( ALinkFileVogtWithWrongSettings, LoadThrowsExceptionIfWrongNdim )
{
	const int Nc = 2;
	const int WrongNDIM = 3;
	LinkFile<PatternStub<float,WrongNDIM,Nc> >* linkfile;
	const int size[4] = {8,4,4,4};
	float* U;
	U = new float[1];

	linkfile = new LinkFileVogt<PatternStub<float,WrongNDIM,Nc>,float>(size);
	linkfile->setFilename( "test_configSU2N4T8SP.vogt" );

	ASSERT_THROW( linkfile->load( U ), LinkFileVogtException );
}

