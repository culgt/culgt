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
#include "LocalLink.h"
#include "GlobalLink.h"
#include "parameterization_types/SUNRealFull.h"
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
	LinkFileVogt( const int size[memoryNdim] ) : LinkFile<MemoryConfigurationPattern>( size )
	{
	}
	LinkFileVogt( const LatticeDimension<memoryNdim> size ) : LinkFile<MemoryConfigurationPattern>( size )
	{
	}

	virtual void loadImplementation() override
	{
		loadHeader();
		verify();
		loadBody();
	}

	void loadHeader()
	{
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&ndim, sizeof(short) );
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&nc, sizeof(short) );
		readSize();
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&sizeOfReal, sizeof(short) );
	}

	LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > getNextLink()
	{
		typedef LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > LocalLink;
		LocalLink link;

		for( int i = 0; i < MemoryConfigurationPattern::PARAMTYPE::SIZE; i++ )
		{
			typename MemoryConfigurationPattern::PARAMTYPE::TYPE value;
			LinkFile<MemoryConfigurationPattern>::file.read( (char*)&value, sizeof(typename MemoryConfigurationPattern::PARAMTYPE::TYPE) );
			link.set( i, value );
		}
		return link;
	}

	void loadBody()
	{
		for( int i = 0; i < this->getLatticeDimension().getSize(); i++ )
		{
			for( int mu = 0; mu < memoryNdim; mu++ )
			{
				typename MemoryConfigurationPattern::SITETYPE site( this->getLatticeDimension() );
				site.setLatticeIndex( i );

				GlobalLink<MemoryConfigurationPattern> dest( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > src;
				src = getNextLink();

				dest = src;
			}
		}
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
			if( this->getLatticeDimension().getDimension(i) != size[i] )
			{
				std::stringstream msg;
				msg << "Wrong lattice size in ";
				msg << i;
				msg << " direction";
				throwException( msg.str(), this->getLatticeDimension().getDimension(i), size[i] );
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

}

#endif /* LINKFILE_H_ */

