/**
 * Very quick and dirty implementation of ILDG gauge config file.
 * Just extract the binary link part and ignore all the rest.
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef LINKFILEILDG_H_
#define LINKFILEILDG_H_

#include <string>
#include <fstream>
#include <sstream>
#include "LinkFile.h"
#include "lattice/LocalLink.h"
#include "lattice/GlobalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "common/culgt_typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "lime_reader.h"
#ifdef __cplusplus
}
#endif

namespace culgt
{

template<typename MemoryConfigurationPattern, typename TFloatFile> class LinkFileILDG: public LinkFile<MemoryConfigurationPattern>
{
private:
	static const lat_dim_t memoryNdim = MemoryConfigurationPattern::SITETYPE::NDIM;

	n_uint64_t headerSize;
	n_uint64_t dataStart;
	n_uint64_t footerStart;
	n_uint64_t footerSize;

	char* header;
	char* footer;
	n_uint64_t allocatedHeaderSize;
	n_uint64_t allocatedFooterSize;

public:
	typedef LinkFile<MemoryConfigurationPattern> super;

	LinkFileILDG(){};
	LinkFileILDG( const int size[memoryNdim], ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
		allocatedHeaderSize = 0;
		allocatedFooterSize = 0;
	}
	LinkFileILDG( const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
		allocatedHeaderSize = 0;
		allocatedFooterSize = 0;
	}

#if __cplusplus == 201103L
	virtual void saveImplementation() override
#else
	void saveImplementation()
#endif
	{
		saveHeader();
		saveBody();
		saveFooter();
	}

#if __cplusplus == 201103L
	virtual void loadImplementation() override
#else
	void loadImplementation()
#endif
	{
		// get the lime record information

		readLimeRecords();

		loadHeader();
		loadBody();
		loadFooter();

//		verify();
	}

	double readDouble()
	{
		int64_t temp;
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&temp, sizeof(int64_t) );
		temp =  __builtin_bswap64( temp );
		double* result = reinterpret_cast<double*>( &temp );
		return *result;
	}

	float readFloat()
	{
		int32_t temp;
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&temp, sizeof(int32_t) );
		temp =  __builtin_bswap32( temp );
		float* result = reinterpret_cast<float*>( &temp );
		return *result;
	}

	void writeDouble( double out )
	{
		int64_t* temp = reinterpret_cast<int64_t*>(&out);
		int64_t result = __builtin_bswap64( *temp );
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&result, sizeof(int64_t) );
	}

	void writeFloat( float out )
	{
		int32_t* temp = reinterpret_cast<int32_t*>(&out);
		int32_t result = __builtin_bswap32( *temp );
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&result, sizeof(int32_t) );
	}

	void readLimeRecords()
	{
		LimeReader* reader;
		FILE *fp;

		fp = fopen( super::getFilename().c_str() , "rb");

		reader = limeCreateReader(fp);

		headerSize = 0;
		footerSize = 0;

		int status;
		char* lime_type;

		while( (status = limeReaderNextRecord(reader)) != LIME_EOF ){

			if( status != LIME_SUCCESS ) {
			  fprintf(stderr, "limeReaderNextRecord returned status = %d\n",
				  status);
			  exit(1);
			}

			n_uint64_t nbytes    = limeReaderBytes(reader);
			lime_type = limeReaderType(reader);

			if( strcmp( lime_type, "ildg-binary-data" ) == 0 )
			{
				headerSize = reader->rec_start;
				footerStart = headerSize + nbytes;
				fseek(fp, 0L, SEEK_END);
				n_uint64_t fileSize = ftell(fp);
				footerSize = fileSize - footerStart;
				break;
			}
		}
		printf("\n");
		printf("headerSize:    %llu\n", (unsigned long long)headerSize);
		printf("footerStart:    %llu\n", (unsigned long long)footerStart);
		printf("footerSize:    %llu\n", (unsigned long long)footerSize);
	}

	void allocateHeader( n_uint64_t size )
	{
		std::cout << "allocate header memory ...";
		if( size > allocatedHeaderSize )
		{
			if( allocatedHeaderSize > 0 ) free( header );
			header = (char*) malloc( size );
			allocatedHeaderSize = size;
		}
		std::cout << "done!" << std::endl;
	}

	void allocateFooter( n_uint64_t size )
	{
		if( size > allocatedFooterSize )
		{
			if( allocatedFooterSize > 0 ) free( footer );
			footer = (char*) malloc( size );
			allocatedFooterSize = size;
		}
	}

	void loadHeader()
	{
		allocateHeader( headerSize );
		LinkFile<MemoryConfigurationPattern>::file.read( header,  headerSize );
		std::cout << "header loaded" << std::endl;
	}

	void loadFooter()
	{
		allocateFooter( footerSize );
		LinkFile<MemoryConfigurationPattern>::file.read( footer,  footerSize );
		std::cout << "footer loaded" << std::endl;
	}

	void saveHeader()
	{
		LinkFile<MemoryConfigurationPattern>::file.write( header, headerSize );
	}

	void saveFooter()
	{
		LinkFile<MemoryConfigurationPattern>::file.write( footer, footerSize );
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > getNextLink()
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> LocalLinkParamType;
		typedef LocalLink<LocalLinkParamType> LocalLink;
		LocalLink link;

		for( int i = 0; i < LocalLinkParamType::SIZE; i++ )
		{
			typename LocalLinkParamType::TYPE value;
			if( super::reinterpretReal == STANDARD )
			{
//				float tempValue;
//				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&tempValue, sizeof(float) );
				value = (typename LocalLinkParamType::TYPE) readFloat();
			}
			else if( super::reinterpretReal == FLOAT )
			{
//				float tempValue;
//				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&tempValue, sizeof(float) );
				value = (typename LocalLinkParamType::TYPE) readFloat();
			}
			else// if( super::reinterpretReal == DOUBLE )
			{
//				double tempValue;
//				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&tempValue, sizeof(double) );
				value = (typename LocalLinkParamType::TYPE) readDouble();
			}
			link.set( i, value );
		}
		return link;
	}

	void writeNextLink( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > link )
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> LocalLinkParamType;

		for( int i = 0; i < LocalLinkParamType::SIZE; i++ )
		{
			if( super::reinterpretReal == STANDARD )
			{
				float value = (float)link.get(i);
				writeFloat( value );
			}
			else if( super::reinterpretReal == FLOAT )
			{
				float value = (float)link.get(i);
				writeFloat( value );
			}
			else// if( super::reinterpretReal == DOUBLE )
			{
				double value = (double)link.get(i);
				writeDouble( value );
			}
		}
	}


	void loadBody()
	{
		SiteCoord<4,NO_SPLIT> siteFile( super::getLatticeDimension() );
		for( int t = 0; t < super::getLatticeDimension().getDimension(0); t++ )
		{
			siteFile[0] = t;
			for( int z = 0; z < super::getLatticeDimension().getDimension(3); z++ )
			{
				siteFile[3] = z;
				for( int y = 0; y < super::getLatticeDimension().getDimension(2); y++ )
				{
					siteFile[2] = y;
					for( int x = 0; x < super::getLatticeDimension().getDimension(1); x++ )
					{
						siteFile[1] = x;
						for( int filemu = 0; filemu < memoryNdim; filemu++ )
						{
							int mu;
							if( filemu == 3 )
							{
								mu = 0;
							}
							else
							{
								mu = filemu+1;
							}
							typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension(), NULL );
							site.setIndexFromNonParitySplitOrder( siteFile.getIndex() );

							GlobalLink<MemoryConfigurationPattern> dest( this->getPointerToU(), site, mu );

							LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > src;
							src = getNextLink();

							dest = src;
						}
					}
				}
			}
		}
	}

	void saveBody()
	{
		SiteCoord<4,NO_SPLIT> siteFile( super::getLatticeDimension() );
		for( int t = 0; t < super::getLatticeDimension().getDimension(0); t++ )
		{
			siteFile[0] = t;
			for( int z = 0; z < super::getLatticeDimension().getDimension(3); z++ )
			{
				siteFile[3] = z;
				for( int y = 0; y < super::getLatticeDimension().getDimension(2); y++ )
				{
					siteFile[2] = y;
					for( int x = 0; x < super::getLatticeDimension().getDimension(1); x++ )
					{
						siteFile[1] = x;
						for( int filemu = 0; filemu < memoryNdim; filemu++ )
						{
							int mu;
							if( filemu == 3 )
							{
								mu = 0;
							}
							else
							{
								mu = filemu+1;
							}


							typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension(), NULL );
							site.setIndexFromNonParitySplitOrder( siteFile.getIndex() );

							GlobalLink<MemoryConfigurationPattern> src( this->getPointerToU(), site, mu );

							LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > dest;

							dest = src;

							writeNextLink( dest );
						}
					}
				}
			}
		}



		for( int i = 0; i < this->getLatticeDimension().getSize(); i++ )
		{
			for( int mu = 0; mu < memoryNdim; mu++ )
			{
				typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension(), NULL );
				site.setIndexFromNonParitySplitOrder( i );

				GlobalLink<MemoryConfigurationPattern> src( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > dest;

				dest = src;

				writeNextLink( dest );
			}
		}
	}

};

}

#endif /* LINKFILE_H_ */

