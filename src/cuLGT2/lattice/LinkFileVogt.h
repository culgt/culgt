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

	void loadBody()
	{
		for( int i = 0; i < latticeSize*memoryNdim; i++ )
		{
			// loadLink
			// assignLink



//			U[MemoryConfigurationPattern::getIndexFromStandardIndex( i )]
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

