/**
 * TODO implement promotion to double / demotion to float
 *
 *
 * ConfigurationFile.hxx
 *
 *  Created on: Apr 16, 2012
 *      Author: vogt
 */

#ifndef CONFIGURATIONFILE_HXX_
#define CONFIGURATIONFILE_HXX_

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../util/datatype/datatypes.h"

template<class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> class LinkFile
{
public:
	LinkFile();
	virtual ~LinkFile();
	bool load( TheSite site, std::string filename, Real *U );
	bool save( std::string filename, Real U );
};

template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::LinkFile()
{
}

template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::~LinkFile()
{
}

/**
 *
 */
template <class Header, class Footer, class FilePattern, class MemoryPattern, class TheSite> bool LinkFile<Header, Footer, FilePattern, MemoryPattern, TheSite>::load( TheSite site, std::string filename, Real *U )
{
	// open file
	std::fstream file;
	file.open( filename.c_str(), std::ios::in | std::ios::binary);

	if( !file )
	{
		util::Logger::log( util::ERROR, "Can't open file");
		util::Logger::log( util::ERROR, filename.c_str() );
		return false;
	}

	// load header

	Header header;
	if( !header.load( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't read header");
	}

	// load config

	// calculate number of reals in configuration:
	int configSize = 1;
	// calculate lattice size
	for( int i = 0; i < site.Ndim; i++ )
	{
		configSize *= site.size[i];
	}

	std::cout << "lattice size: " << configSize << std::endl;

	// * number of directions mu * Nc^2 * 2 (complex number)
	configSize *= FilePattern::Ndim * FilePattern::Nc * FilePattern::Nc * 2;

	std::cout << "Ndim: " << FilePattern::Ndim << std::endl;

	std::cout << "total size: " << configSize << std::endl;

	Real temp;
	for( int i = 0; i < configSize; i++ )
	{
		// load number
		file.read( (char*)&temp, sizeof(Real) );
//		std::cout << temp << std::endl;
//		if( i < 18 ) std::cout << temp << std::endl;
		int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
		U[index] = temp;
//		break;
	}


	// load footer
//	if( !footer.load( file ) )
//	{
//		util::Logger::log( util::ERROR, "Can't read footer");
//	}


	file.close();
	return true;
}


#endif /* CONFIGURATIONFILE_HXX_ */
