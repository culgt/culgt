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

namespace culgt
{


class IOException: public std::exception
{
public:
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

	virtual void load( typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U ) final
	{
		this->U = U;
		this->loadImplementation();
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

protected:
	virtual void loadImplementation(){};

private:
	std::string filename;
	std::fstream file;
	bool fileIsOpen;
	bool filenameIsSet;
	typename MemoryConfigurationPattern::PARAMTYPE::TYPE* U;
};

}

#endif /* LINKFILE_H_ */


