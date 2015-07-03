/**
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef LINKFILENERSC_H_
#define LINKFILENERSC_H_

#include <string>
#include <fstream>
#include <sstream>
#include "LinkFile.h"
#include "lattice/LocalLink.h"
#include "lattice/GlobalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/SU3Real12.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU3_Real12_Real18.h"
#include "common/culgt_typedefs.h"
#include "common/culgt_compat.h"
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>

using std::string;
using std::istringstream;

namespace culgt
{




enum NerscDatatype { NERSC_4D_SU3_GAUGE_3x3, NERSC_4D_SU3_GAUGE };

template<typename MemoryConfigurationPattern> class LinkFileNERSC: public LinkFile<MemoryConfigurationPattern>
{
private:
	typedef typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE REALTYPE;
	static const lat_dim_t memoryNdim = MemoryConfigurationPattern::SITETYPE::NDIM;

	std::map<string,string> headerData;

	bool bigendian;
	bool isDouble;
	int latticeSizeInFile[4];
	NerscDatatype datatype;

	bool fullHeaderIsAvailable;

public:
	typedef LinkFile<MemoryConfigurationPattern> super;

	LinkFileNERSC(): datatype(NERSC_4D_SU3_GAUGE), fullHeaderIsAvailable(false)
	{};
	LinkFileNERSC( const int size[memoryNdim], ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret ), datatype(NERSC_4D_SU3_GAUGE), fullHeaderIsAvailable(false)
	{
	}
	LinkFileNERSC( const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret ), datatype(NERSC_4D_SU3_GAUGE), fullHeaderIsAvailable(false)
	{
	}

	virtual void saveImplementation() CULGT_OVERRIDE
	{
		saveHeader();
		saveBody();
	}

	virtual void loadImplementation() CULGT_OVERRIDE
	{
		loadHeader();
		verify();
		loadBody();
	}

	std::vector<string> readHeaderLines()
	{
		std::vector<string> headerLines;
		string line;
		while( std::getline( super::file, line ) )
		{
			headerLines.push_back( line );
			if( line.compare( "END_HEADER" ) == 0 ) break;
		}
		return headerLines;
	}

	bool isNotBeginOrEnd( string key )
	{
		return key.compare("BEGIN_HEADER") != 0 && key.compare("END_HEADER") != 0;
	}

	void parseFloatType()
	{
		string floatingPoint = headerData["FLOATING_POINT"];
		if( floatingPoint.compare( "IEEE32" ) == 0 )
		{
			isDouble = false;
			bigendian = false;
		}
		else if( floatingPoint.compare( "IEEE32BIG" ) == 0 )
		{
			isDouble = false;
			bigendian = true;
		}
		else if( floatingPoint.compare( "IEEE64" ) == 0 )
		{
			isDouble = true;
			bigendian = false;
		}
		else if( floatingPoint.compare( "IEEE64BIG" ) == 0 )
		{
			isDouble = true;
			bigendian = true;

		}
		else
		{
			throw LinkFileException( "unknown floating point format" );
		}
	}

	string getFloatingPointString()
	{
		if( getSizeOfReal() == 8 )
		{
			if( bigendian )
				return "IEEE64BIG";
			else
				return "IEEE64";
		}
		else
		{
			if( bigendian )
				return "IEEE32BIG";
			else
				return "IEEE32";
		}
	}

	void parseLatticeSize()
	{
		istringstream( headerData["DIMENSION_4"] ) >> latticeSizeInFile[0];
		istringstream( headerData["DIMENSION_1"] ) >> latticeSizeInFile[1];
		istringstream( headerData["DIMENSION_2"] ) >> latticeSizeInFile[2];
		istringstream( headerData["DIMENSION_3"] ) >> latticeSizeInFile[3];
		for( int i = 0; i < 4; i++ )
			std::cout << latticeSizeInFile[i] << " ";
		std::cout << std::endl;
	}

	void parseNerscDatatype()
	{
		if( headerData["DATATYPE"].compare("4D_SU3_GAUGE_3x3") == 0 )
			datatype = NERSC_4D_SU3_GAUGE_3x3;
		else if( headerData["DATATYPE"].compare("4D_SU3_GAUGE") == 0 )
			datatype = NERSC_4D_SU3_GAUGE;
		else throw LinkFileException( "unknown NERSC datatype" );
	}

	void parseHeader( std::vector<string>& headerLines )
	{
		for( std::vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); ++it )
		{
			istringstream stream( *it );
			string key;
			std::getline( stream, key, '=' );
			boost::trim( key );

			if( isNotBeginOrEnd( key ) )
			{
				string value;
				std::getline( stream, value );
				boost::trim( value );

				headerData[key] = value;
			}
		}

		parseFloatType();
		parseLatticeSize();
		parseNerscDatatype();
	}

	void loadHeader()
	{
		std::vector<string> headerLines = readHeaderLines();
		parseHeader( headerLines );

		fullHeaderIsAvailable = true;

		for( std::map<string,string>::iterator it = headerData.begin(); it != headerData.end(); ++it )
		{
			std::cout << it->first << " = " << it->second << std::endl;
		}
	}

	void makeHeaderData()
	{
		headerData.clear();
		headerData["HDR_VERSION"] = "1.0";
		if( datatype == NERSC_4D_SU3_GAUGE )
		{
			headerData["DATATYPE"] = "4D_SU3_GAUGE";
		}
		else
		{
			headerData["DATATYPE"] = "4D_SU3_GAUGE_3x3";
		}
		headerData["DIMENSION_1"] = to_string( this->getLatticeDimension().getDimension(1) );
		headerData["DIMENSION_2"] = to_string( this->getLatticeDimension().getDimension(2) );
		headerData["DIMENSION_3"] = to_string( this->getLatticeDimension().getDimension(3) );
		headerData["DIMENSION_4"] = to_string( this->getLatticeDimension().getDimension(0) );


		headerData["FLOATING_POINT"] = getFloatingPointString();

		headerData["CHECKSUM"] = "0";
		headerData["LINK_TRACE"] = "0";
		headerData["PLAQUETTE"] = "0";
		headerData["CREATOR"] = "n/a";
		headerData["CREATOR_HARDWARE"] = "n/a";
		headerData["ARCHIVE_DATE"] = "Thu Jan 01 00:00:00 1970";
		headerData["ENSEMBLE_LABEL"] = "n/a";
		headerData["ENSEMBLE_ID"] = "000";
		headerData["BETA"] = "0.0";
	}

	void saveHeader()
	{
		if( !fullHeaderIsAvailable )
		{
			makeHeaderData();
		}

		super::file << "BEGIN_HEADER" << std::endl;
		for( std::map<string,string>::iterator it = headerData.begin(); it != headerData.end(); ++it )
		{
			super::file << it->first << " = " << it->second << std::endl;
		}
		super::file << "END_HEADER" << std::endl;
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > getNextLink()
	{
		if( datatype == NERSC_4D_SU3_GAUGE_3x3 ) return getNextLink3x3();
		else return getNextLink3x2();
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > getNextLink3x3()
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> LocalLinkParamType;
		typedef LocalLink<LocalLinkParamType> LocalLink;
		LocalLink link;

		for( int i = 0; i < 18; i++ )
		{
			typename LocalLinkParamType::TYPE value;

			if( isDouble )
				value = (typename LocalLinkParamType::TYPE) super::readDouble( bigendian );
			else
				value = (typename LocalLinkParamType::TYPE) super::readFloat( bigendian );

			link.set( i, value );
		}
		return link;
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > getNextLink3x2()
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> FullLinkParamType;
		typedef LocalLink<FullLinkParamType> FullLink;

		typedef SU3Real12<REALTYPE> LocalLinkParamType;
		typedef LocalLink<LocalLinkParamType> LocalLink;
		LocalLink link;

		for( int i = 0; i < 12; i++ )
		{
			typename LocalLinkParamType::TYPE value;

			if( isDouble )
				value = (typename LocalLinkParamType::TYPE) super::readDouble( bigendian );
			else
				value = (typename LocalLinkParamType::TYPE) super::readFloat( bigendian );

			link.set( i, value );
		}

		FullLink fullLink;
		fullLink = link;
		return fullLink;
	}

	void writeNextLink( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > link )
	{
		if( datatype == NERSC_4D_SU3_GAUGE_3x3 ) writeNextLink3x3( link );
			else writeNextLink3x2( link );
	}

	void writeNextLink3x3( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > link )
	{
		writeNextLink( link, 3*3*2 );
	}

	void writeNextLink3x2( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > link )
	{
		writeNextLink( link, 3*2*2 );
	}

	void writeNextLink( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > link, int parametersToWrite )
	{
		for( int i = 0; i < parametersToWrite; i++ )
		{
			if( getSizeOfReal() == 8 )
			{
				double value = (double)link.get(i);
				super::writeDouble( value, bigendian );
			}
			else
			{
				float value = (float)link.get(i);
				super::writeFloat( value, bigendian );
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

							LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > src;
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

							LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > dest;

							dest = src;

							writeNextLink( dest );
						}
					}
				}
			}
		}
	}

	int getSizeOfReal()
	{
		if( super::reinterpretReal == FLOAT ) return sizeof( float );
		else if( super::reinterpretReal == DOUBLE ) return sizeof( double );
		else return sizeof( typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE );
	}

	void verify()
	{
		int mySizeOfReal = getSizeOfReal();

		int fileSizeOfReal;
		if( isDouble ) fileSizeOfReal = 8;
		else fileSizeOfReal = 4;

		if( mySizeOfReal != fileSizeOfReal )
		{
			super::throwException( "Wrong size of real", mySizeOfReal, fileSizeOfReal );
		}

		for( int i = 0; i < 4; i++ )
		{
			if( this->getLatticeDimension().getDimension(i) != latticeSizeInFile[i] )
			{
				std::stringstream msg;
				msg << "Wrong lattice size in ";
				msg << i;
				msg << " direction";
				super::throwException( msg.str(), this->getLatticeDimension().getDimension(i), latticeSizeInFile[i] );
			}
		}
	}
};

}

#endif /* LINKFILE_H_ */

