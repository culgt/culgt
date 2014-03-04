/**
 * LinkFile.h
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef LINKFILEVOGT_H_
#define LINKFILEVOGT_H_

#include <string>
#include <fstream>
#include <sstream>
#include "LinkFile.h"
#include "../common/culgt_typedefs.h"

namespace culgt
{


class LinkFileVogtException: public IOException
{
public:
	LinkFileVogtException( std::string msg ) : IOException(msg){};
};




template<typename MemoryConfigurationPattern, typename TFloatFile> class LinkFileVogt: public LinkFile<MemoryConfigurationPattern>
{
private:
	short ndim;
	short nc;
	static const lat_dim_t memoryNdim = MemoryConfigurationPattern::SITETYPE::Ndim;
	short size[memoryNdim];
	short memorySize[memoryNdim];
	short sizeOfReal;

	void readSize()
	{
		for( int i = 0; i < memoryNdim; i++ )
		{
			LinkFile<MemoryConfigurationPattern>::file.read( (char*)&size[i], sizeof(short) );
		}
	}

	void throwException( std::string msg, short expected, short inFile )
	{
		std::stringstream message;
		message << msg;
		message << ": Expected ";
		message << expected;
		message << " in file, got ";
		message << inFile;
		throw LinkFileVogtException( message.str() );
	}

public:
	LinkFileVogt(){};
	LinkFileVogt( const int size[memoryNdim] )
	{
		for( int i = 0; i < memoryNdim; i++ )
		{
			this->memorySize[i] = size[i];
		}
	}

	virtual void loadImplementation() override;
	void loadHeader()
	{
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&ndim, sizeof(short) );
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&nc, sizeof(short) );
		readSize();
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&sizeOfReal, sizeof(short) );
	}

	void verify()
	{
		if( memoryNdim != ndim )
		{
			throwException( "Wrong lattice dimension", (short)memoryNdim, ndim );
		}

		if( MemoryConfigurationPattern::PARAMTYPE::NC != nc )
		{
			throwException( "Wrong gauge group", MemoryConfigurationPattern::PARAMTYPE::NC, nc );
		}

		if( sizeof( typename MemoryConfigurationPattern::PARAMTYPE::TYPE ) != sizeOfReal )
		{
			throwException( "Wrong size of real", sizeof( typename MemoryConfigurationPattern::PARAMTYPE::TYPE ), sizeOfReal );
		}

		for( int i = 0; i < memoryNdim; i++ )
		{
			if( memorySize[i] != size[i] )
			{
				std::stringstream msg;
				msg << "Wrong lattice size in ";
				msg << i;
				msg << " direction";
				throwException( msg.str(), memorySize[i], size[i] );
			}
		}
	}

	void loadData()
	{

	}

	short getNdim() const
	{
		return ndim;
	}

	short getNc() const
	{
		return nc;
	}

	short getSize( int i ) const
	{
		return size[i];
	}

	short getSizeOfReal() const
	{
		return sizeOfReal;
	}

};



template<typename MemoryConfigurationPattern, typename TFloatFile > void LinkFileVogt<MemoryConfigurationPattern, TFloatFile>::loadImplementation()
{
	loadHeader();
	verify();
}


}

#endif /* LINKFILE_H_ */


#include "gmock/gmock.h"
#include "testhelper_pattern_stub.h"

using namespace culgt;
using namespace ::testing;

/**
 * The sample configuration is 8x4^3 lattice SU(2) lattice in Single Precision with first entry -1.966492e-01
 */
class ALinkFileVogtWithSampleConfiguration: public Test
{
public:
	LinkFileVogt<PatternStub<float>,float> linkfile;

	void SetUp()
	{
		linkfile.setFilename("../test_configSU2N4T8SP.vogt");
		linkfile.openFile();
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

TEST_F( ALinkFileVogtWithSampleConfiguration, FirstEntryAfter )

template<typename LinkFileType> void readSampleFileHeader( LinkFileType& linkfile )
{
	linkfile.setFilename("../test_configSU2N4T8SP.vogt");
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
	linkfile->setFilename( "../test_configSU2N4T8SP.vogt" );

	ASSERT_THROW( linkfile->load( U ), LinkFileVogtException );
}

