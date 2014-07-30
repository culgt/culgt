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

	void writeSize()
	{
		for( int i = 0; i < memoryNdim; i++ )
		{
			short mySize = this->getLatticeDimension().getDimension(i);
			LinkFile<MemoryConfigurationPattern>::file.write( (char*)&mySize, sizeof(short) );
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
	typedef LinkFile<MemoryConfigurationPattern> super;

	LinkFileVogt(){};
	LinkFileVogt( const int size[memoryNdim], ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
	}
	LinkFileVogt( const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret = STANDARD ) : LinkFile<MemoryConfigurationPattern>( size, reinterpret )
	{
	}

#if __cplusplus == 201103L
	virtual void saveImplementation() override
#else
	void saveImplementation()
#endif
	{
		saveHeader();
		saveBody();
	}

#if __cplusplus == 201103L
	virtual void loadImplementation() override
#else
	void loadImplementation()
#endif
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

	void saveHeader()
	{
		short myMemoryNdim = memoryNdim;
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&myMemoryNdim, sizeof(short) );
		short myNc= MemoryConfigurationPattern::PARAMTYPE::NC;
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&myNc, sizeof(short) );
		writeSize();
		short mySizeOfReal;
		if( super::reinterpretReal == STANDARD) mySizeOfReal = sizeof( typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE );
		else if( super::reinterpretReal == FLOAT ) mySizeOfReal = sizeof( float );
		else if( super::reinterpretReal == DOUBLE ) mySizeOfReal = sizeof( double );
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&mySizeOfReal, sizeof(short) );
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

	void writeNextLink( LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > link )
	{
		typedef SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> LocalLinkParamType;

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
				site.setLatticeIndexFromNonParitySplitOrder( i );

				GlobalLink<MemoryConfigurationPattern> dest( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > src;
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
				site.setLatticeIndexFromNonParitySplitOrder( i );

				GlobalLink<MemoryConfigurationPattern> src( this->getPointerToU(), site, mu );

				LocalLink<SUNRealFull<MemoryConfigurationPattern::PARAMTYPE::NC,TFloatFile> > dest;

				dest = src;

				writeNextLink( dest );
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


		short mySizeOfReal;
		if( super::reinterpretReal == STANDARD) mySizeOfReal = sizeof( typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE );
		else if( super::reinterpretReal == FLOAT ) mySizeOfReal = sizeof( float );
		else if( super::reinterpretReal == DOUBLE ) mySizeOfReal = sizeof( double );

		if( mySizeOfReal != sizeOfReal )
		{
			throwException( "Wrong size of real", sizeof( typename MemoryConfigurationPattern::PARAMTYPE::REALTYPE ), sizeOfReal );
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

