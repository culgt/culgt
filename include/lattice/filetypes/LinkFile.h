/**
 * LinkFile.h
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef LINKFILE_H_
#define LINKFILE_H_

#include <string>
#include <fstream>
#include "lattice/LatticeDimension.h"
#include "common/culgt_typedefs.h"
#include "common/culgt_compat.h"
#include "filetype_typedefs.h"

namespace culgt
{

class IOException: public std::exception
{
public:
	~IOException() throw() {};
	IOException( std::string msg ) : msg(msg){};
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
protected:
	std::string msg;
};

class LinkFileException: public IOException
{
public:
	LinkFileException( std::string msg ) : IOException(msg){};
};

class FileNotFoundException: IOException
{
public:
	FileNotFoundException( std::string filename ) : IOException( "File " + filename+ " not found!" ) {};
};



template<typename MemoryConfigurationPattern> class LinkFile
{
public:
	virtual ~LinkFile() {};


	LinkFile() : reinterpretReal(STANDARD), fileIsOpen( false ), filenameIsSet( false )
	{
	};

	LinkFile( ReinterpretReal reinterpret ) : reinterpretReal(reinterpret), fileIsOpen( false ), filenameIsSet( false )
	{
	};

	LinkFile( const int size[MemoryConfigurationPattern::SITETYPE::NDIM], ReinterpretReal reinterpret ) :reinterpretReal(reinterpret), fileIsOpen( false ), filenameIsSet( false ), latticeDimension( size )
	{
	};

	LinkFile( const LatticeDimension<MemoryConfigurationPattern::SITETYPE::NDIM> size, ReinterpretReal reinterpret ) : reinterpretReal(reinterpret),fileIsOpen( false ), filenameIsSet( false ), latticeDimension( size )
	{
	};

	double readDouble( bool convertBigEndian = true )
	{
		int64_t temp;
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&temp, sizeof(int64_t) );
		if( convertBigEndian ) temp =  __builtin_bswap64( temp );
		double* result = reinterpret_cast<double*>( &temp );
		return *result;
	}

	float readFloat( bool convertBigEndian = true )
	{
		int32_t temp;
		LinkFile<MemoryConfigurationPattern>::file.read( (char*)&temp, sizeof(int32_t) );
		if( convertBigEndian ) temp =  __builtin_bswap32( temp );
		float* result = reinterpret_cast<float*>( &temp );
		return *result;
	}

	void writeDouble( double out, bool convertBigEndian = true  )
	{
		int64_t* temp = reinterpret_cast<int64_t*>(&out);
		int64_t result;
		if( convertBigEndian )
			result = __builtin_bswap64( *temp );
		else
			result = *temp;
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&result, sizeof(int64_t) );
	}

	void writeFloat( float out, bool convertBigEndian = true  )
	{
		int32_t* temp = reinterpret_cast<int32_t*>(&out);
		int32_t result;
		if( convertBigEndian )
			result = __builtin_bswap32( *temp );
		else
			result = *temp;
		LinkFile<MemoryConfigurationPattern>::file.write( (char*)&result, sizeof(int32_t) );
	}

	virtual void load( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U ) CULGT_FINAL
	{
		this->U = U;
		this->openFile();
		this->loadImplementation();
		this->closeFile();
	};

	virtual void save( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U ) CULGT_FINAL
	{
		this->U = U;
		this->openFileWrite();
		this->saveImplementation();
		this->closeFile();
	};

	void openFileWrite()
	{
		if( !filenameIsSet )
		{
			throw IOException( "Filename not set." );
		}

		file.open( filename.c_str(), std::ios::out | std::ios::binary);

		if( file.good() )
		{
			fileIsOpen = true;
		}
		else
		{
			throw FileNotFoundException( filename );
		}
	};

	void openFile()
	{
		if( !filenameIsSet )
		{
			throw IOException( "Filename not set." );
		}

		file.open( filename.c_str(), std::ios::in | std::ios::binary);

		if( file.good() )
		{
			fileIsOpen = true;
		}
		else
		{
			throw FileNotFoundException( filename );
		}
	};
	void closeFile()
	{
		if( fileIsOpen )
		{
			file.close();
		}
	}
	const std::string& getFilename() const
	{
		return filename;
	}

	void setFilename(const std::string& filename)
	{
		this->filename = filename;
		filenameIsSet = true;
	}

	void setPointerToU( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U )
	{
		this->U = U;
	}

	typename MemoryConfigurationPattern::PARAMTYPE::TYPE* getPointerToU() const
	{
		return U;
	}

	LatticeDimension<MemoryConfigurationPattern::SITETYPE::NDIM> getLatticeDimension() const
	{
		return latticeDimension;
	}


protected:
	std::fstream file;
	ReinterpretReal reinterpretReal;
	virtual void loadImplementation() = 0;
	virtual void saveImplementation() = 0;
	virtual std::string getPreferredExtension() = 0;

	void throwException( std::string msg, int expected, int inFile )
	{
		std::stringstream message;
		message << msg;
		message << ": Expected ";
		message << expected;
		message << " in file, got ";
		message << inFile;
		throw LinkFileException( message.str() );
	}

private:
	bool fileIsOpen;
	bool filenameIsSet;
	LatticeDimension<MemoryConfigurationPattern::SITETYPE::NDIM> latticeDimension;
	std::string filename;
	typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U;
};

}

#endif /* LINKFILE_H_ */


