/**
 * Class to read configurations in a "generic form" of the mdp-format,
 * i.e., the same data ordering as in the mdp-format is expected
 * but arbitrary headers are allowed (will simply be copied to the output file).
 * No footer allowed!
 *
 * TODO give an explict definition of th format here!
 */

#ifndef LINKFILEHEADERONLY_H_
#define LINKFILEHEADERONLY_H_

#include <string>
#include <fstream>
#include <sstream>
#include "LinkFile.h"
#include "lattice/LocalLink.h"
#include "lattice/GlobalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "common/culgt_typedefs.h"
#include "common/culgt_compat.h"
#include <iosfwd>

using std::ios;

namespace culgt
{

template<typename MemoryConfigurationPattern> class LinkFileHeaderOnly: public LinkFile<MemoryConfigurationPattern>
{
private:
	static const lat_dim_t memoryNdim = MemoryConfigurationPattern::SITETYPE::NDIM;
	typedef typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE REALTYPE;

	long arraySize;
	long filelength;
	int offset;
	char *header;
public:
	typedef LinkFile<MemoryConfigurationPattern> super;
	typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> LocalLinkParamType;

	LinkFileHeaderOnly(){};
	LinkFileHeaderOnly( const int size[memoryNdim], ReinterpretReal reinterpret ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
		init();
	}
	LinkFileHeaderOnly( const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
		init();
	}

	void init()
	{
		int realsize;
		if( super::reinterpretReal == FLOAT )
		{
			realsize = sizeof(float);
		}
		else if( super::reinterpretReal == DOUBLE )
		{
			realsize = sizeof(double);
		}
		else // STANDARD
		{
			realsize = sizeof(typename LocalLinkParamType::TYPE);
		}

		arraySize = super::getLatticeDimension().getSize()*memoryNdim*LocalLinkParamType::SIZE*realsize;
	}

	virtual void saveImplementation() CULGT_OVERRIDE
	{
		saveHeader();
		saveBody();
	}

	virtual void loadImplementation() CULGT_OVERRIDE
	{
		loadHeader();
		loadBody();
	}

	virtual std::string getPreferredExtension() CULGT_OVERRIDE
	{
		return "";
	}

	void loadHeader()
	{
		//what's the offset (the header length)?
		LinkFile<MemoryConfigurationPattern>::file.seekg(0, ios::end);
	  filelength = LinkFile<MemoryConfigurationPattern>::file.tellg();
	  LinkFile<MemoryConfigurationPattern>::file.seekg(0, ios::beg);
		offset = filelength - arraySize;
		std::cout << "\nfilelength = " << filelength;
		std::cout << "\narraysize = " << arraySize;
		std::cout << "\noffset by header = " << offset << std::endl;

		header = (char*) malloc(offset);

	  //copy header
		LinkFile<MemoryConfigurationPattern>::file.read(header,offset);

		std::cout << header << std::endl;
	}

	void saveHeader()
	{
		LinkFile<MemoryConfigurationPattern>::file.write(header,offset);
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,REALTYPE> > getNextLink()
	{
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
};

}

#endif /* LINKFILE_H_ */

