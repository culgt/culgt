/**
 */

#ifndef LINKFILEMDP_H_
#define LINKFILEMDP_H_

#include <string>
#include <fstream>
#include <sstream>
#include "LinkFile.h"
#include "lattice/LocalLink.h"
#include "lattice/GlobalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "common/culgt_typedefs.h"
#include "common/culgt_compat.h"
#include <boost/algorithm/string.hpp>

using boost::algorithm::trim;

namespace culgt
{

template<typename MemoryConfigurationPattern> class LinkFileMDP: public LinkFile<MemoryConfigurationPattern>
{
private:
	typedef typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE REALTYPE;
	static const lat_dim_t memoryNdim = MemoryConfigurationPattern::SITETYPE::NDIM;

	static const int MAX_DIM = 10;
	static const std::string MDP_FILETYPE_STRING;
	static const unsigned int MDP_ENDIANESS;

	int ndim;
	int size[MAX_DIM];
	int siteSize;

	void readSize()
	{
		for( int i = 0; i < MAX_DIM; i++ )
		{
			super::file.read( (char*)&size[i], sizeof(int) );
		}
	}

	void writeSize()
	{
		for( int i = 0; i < memoryNdim; i++ )
		{
			int mySize = this->getLatticeDimension().getDimension(i);
			LinkFile<MemoryConfigurationPattern>::file.write( (char*)&mySize, sizeof(int) );
		}
		for( int i = memoryNdim; i < MAX_DIM; i++ )
		{
			int mySize = 0;
			LinkFile<MemoryConfigurationPattern>::file.write( (char*)&mySize, sizeof(int) );
		}
	}

	int calculateSiteSize()
	{
		int mySizeOfReal;
		if( super::reinterpretReal == FLOAT ) mySizeOfReal = sizeof( float );
		else if( super::reinterpretReal == DOUBLE ) mySizeOfReal = sizeof( double );
		else mySizeOfReal = sizeof( typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE );

		return MemoryConfigurationPattern::PARAMTYPE::NC*MemoryConfigurationPattern::PARAMTYPE::NC*2*memoryNdim*mySizeOfReal;
	}

public:
	typedef LinkFile<MemoryConfigurationPattern> super;

	LinkFileMDP(){};
	LinkFileMDP( const int size[memoryNdim], ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
	}
	LinkFileMDP( const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
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

	virtual std::string getPreferredExtension() CULGT_OVERRIDE
	{
		return ".mdp";
	}

	void loadHeader()
	{
		char temp[60];

		super::file.read( temp, 60*sizeof(char) );
		std::string filetypeString( temp );
		trim( filetypeString );
		if( filetypeString.compare( MDP_FILETYPE_STRING ) != 0 ) throw LinkFileException( "File header does not contain MDP header string" );

		super::file.read( temp, 60*sizeof(char) ); // skip filename
		super::file.read( temp, 60*sizeof(char) ); // skip date

		unsigned int magicnumber;
		super::file.read( (char*)&magicnumber, sizeof(unsigned int) );
		if( magicnumber != MDP_ENDIANESS ) throw LinkFileException( "MDP endianess check failed" );

		super::file.read( (char*)&ndim, sizeof(int) );

		readSize();

		super::file.read( (char*)&siteSize, sizeof(int) );

		int latticeVolume;
		super::file.read( (char*)&latticeVolume, sizeof(int) ); // skip latticeVolume

	}

	void saveHeader()
	{
		char temp[60];
		sprintf( temp, "%s", MDP_FILETYPE_STRING.c_str());
		super::file.write( temp, 60*sizeof(char) );

		sprintf( temp, "n/a" );
		super::file.write( temp, 60*sizeof(char) ); // write empty filename
		super::file.write( temp, 60*sizeof(char) ); // write empty date

		super::file.write( (char*)&MDP_ENDIANESS, sizeof( unsigned int) );

		super::file.write( (char*)&memoryNdim, sizeof( int) );

		writeSize();

		int siteSize = calculateSiteSize();
		super::file.write( (char*)&siteSize, sizeof( int) );

		int latticeVolume = this->getLatticeDimension().getSize();
		super::file.write( (char*)&latticeVolume, sizeof( int) );
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > getNextLink()
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> LocalLinkParamType;
		typedef LocalLink<LocalLinkParamType> LocalLink;
		LocalLink link;

		for( int i = 0; i < LocalLinkParamType::SIZE; i++ )
		{
			typename LocalLinkParamType::TYPE value;
			if( super::reinterpretReal == STANDARD )
			{
				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&value, sizeof(typename LocalLinkParamType::TYPE) );
			}
			else if( super::reinterpretReal == FLOAT )
			{
				float tempValue;
				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&tempValue, sizeof(float) );
				value = (typename LocalLinkParamType::TYPE) tempValue;
			}
			else if( super::reinterpretReal == DOUBLE )
			{
				double tempValue;
				LinkFile<MemoryConfigurationPattern>::file.read( (char*)&tempValue, sizeof(double) );
				value = (typename LocalLinkParamType::TYPE) tempValue;
			}
			link.set( i, value );
		}
		return link;
	}

	void writeNextLink( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > link )
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> LocalLinkParamType;

		for( int i = 0; i < LocalLinkParamType::SIZE; i++ )
		{
			if( super::reinterpretReal == STANDARD )
			{
				typename LocalLinkParamType::TYPE value = link.get( i );
				LinkFile<MemoryConfigurationPattern>::file.write( (char*)&value, sizeof(typename LocalLinkParamType::TYPE) );
			}
			else if( super::reinterpretReal == FLOAT )
			{
				float value = (float)link.get(i);
				LinkFile<MemoryConfigurationPattern>::file.write( (char*)&value, sizeof(float) );
			}
			else if( super::reinterpretReal == DOUBLE )
			{
				double value = (double)link.get(i);
				LinkFile<MemoryConfigurationPattern>::file.write( (char*)&value, sizeof(double) );
			}
		}
	}


	void loadBody()
	{
		for( int i = 0; i < this->getLatticeDimension().getSize(); i++ )
		{
			for( int mu = 0; mu < memoryNdim; mu++ )
			{
				typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension(), NULL );
				site.setIndexFromNonParitySplitOrder( i );

				GlobalLink<MemoryConfigurationPattern> dest( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > src;
				src = getNextLink();

				dest = src;
			}
		}
	}

	void saveBody()
	{
		for( int i = 0; i < this->getLatticeDimension().getSize(); i++ )
		{
			for( int mu = 0; mu < memoryNdim; mu++ )
			{
				typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension(), NULL );
				site.setIndexFromNonParitySplitOrder( i );

				GlobalLink<MemoryConfigurationPattern> src( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > dest;

				dest = src;

				writeNextLink( dest );
			}
		}
	}

	void verify()
	{
		if( memoryNdim != ndim )
		{
			super::throwException( "Wrong lattice dimension", (short)memoryNdim, ndim );
		}

		for( int i = 0; i < memoryNdim; i++ )
		{
			if( this->getLatticeDimension().getDimension(i) != size[i] )
			{
				std::stringstream msg;
				msg << "Wrong lattice size in ";
				msg << i;
				msg << " direction";
				super::throwException( msg.str(), this->getLatticeDimension().getDimension(i), size[i] );
			}
		}

		if( calculateSiteSize() != siteSize )
		{
			super::throwException( "Wrong site size", calculateSiteSize(), siteSize );
		}

	}
};

template<typename Pattern> const std::string LinkFileMDP<Pattern>::MDP_FILETYPE_STRING = "File Type: MDP FIELD";
template<typename Pattern> const unsigned int LinkFileMDP<Pattern>::MDP_ENDIANESS = 0x87654321;
template<typename Pattern> const lat_dim_t LinkFileMDP<Pattern>::memoryNdim;

}

#endif /* LINKFILE_H_ */

