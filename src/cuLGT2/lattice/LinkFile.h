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
#include "LatticeDimension.h"
#include "../common/culgt_typedefs.h"

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

class FileNotFoundException: IOException
{
public:
	FileNotFoundException( std::string filename ) : IOException( "File " + filename+ " not found!" ) {};
};



template<typename MemoryConfigurationPattern> class LinkFile
{
public:
	virtual ~LinkFile() {};


	LinkFile() : fileIsOpen( false ), filenameIsSet( false )
	{
	};

//	LinkFile( const LinkFile<MemoryConfigurationPattern>& copy ) : fileIsOpen( false ), filenameIsSet( false ), latticeDimension(copy.getLatticeDimension() )
//	{
//		U = copy.getPointerToU();
//	};

	LinkFile( const int size[MemoryConfigurationPattern::SITETYPE::Ndim] ) : fileIsOpen( false ), filenameIsSet( false ), latticeDimension( size )
	{
	};

	LinkFile( const LatticeDimension<MemoryConfigurationPattern::SITETYPE::Ndim> size ) : fileIsOpen( false ), filenameIsSet( false ), latticeDimension( size )
	{
	};

#if __cplusplus == 201103L
	virtual void load( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U ) final
#else
	void load( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U )
#endif
	{
		this->U = U;
		this->openFile();
		this->loadImplementation();
		this->closeFile();
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

	LatticeDimension<MemoryConfigurationPattern::SITETYPE::Ndim> getLatticeDimension() const
	{
		return latticeDimension;
	}


protected:
	std::fstream file;
	virtual void loadImplementation(){};

private:
	bool fileIsOpen;
	bool filenameIsSet;
	LatticeDimension<MemoryConfigurationPattern::SITETYPE::Ndim> latticeDimension;
	std::string filename;
	typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U;
};

}

#endif /* LINKFILE_H_ */


