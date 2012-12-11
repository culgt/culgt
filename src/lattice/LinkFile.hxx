/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 *
 * Class for loading and storing configurations.
 * TODO I named this class LinkFile because ConfigurationFile may be misleading but I don't like the naming.
 *      Suggestions for a more precise name are very welcome.
 *
 * This class serves as a middleman between different array patterns on file and memory side.
 * The pattern for saving configurations is usually such that matrices are stored in contiguous regions the memory pattern
 * for CUDA is completely different. This template class accepts as arguments the FilePattern, i.e. the array pattern on file side,
 * the MemoryPattern.
 * To allow a high flexibility you can define FileType classes that specify header and/or footer that are called before/after the read of the configuration.
 *
 * The intension of this class is to cover all possible File- and MemoryPatterns. This causes some computation overhead which can
 * be avoided by spezializations of this class.
 * Check the getUnique()/setUnique() function calls in load() and save(): Each Pattern needs to be able to provide a unique index
 * which does not depend on the pattern. Calculating back and forth this index maybe costly.
 * TODO Check how much it costs actually.
 *
 * TODO:
 *  - Is this class really flexible? Do we have to move the loading of the configuration to an extra class?
 *  - Implement a switch for promotion to double / demotion to float
 *  - Think about partially loading when we want to use the reconstruction technique for SU3.
 *
 */

#ifndef LINKFILE_HXX_
#define LINKFILE_HXX_

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "datatype/datatypes.h"
#include "../util/log/Logger.hxx"
#include "filetypes/filetype_typedefs.h"


// TODO IMPORTANT: replace TheSite by an array that contains the lattice extents (we must not use TheSite, but only its size, because the types of Site may be different in File- and MemoryPattern)

template<class FileType, class FilePattern, class MemoryPattern, class TheSite> class LinkFile
{
public:
	LinkFile( ReinterpretReal reinterpret = STANDARD );
	virtual ~LinkFile();
	bool load( TheSite site, std::string filename, Real *U );
	bool save( TheSite site, std::string filename, Real *U );
	FileType filetype;
private:
	ReinterpretReal reinterpret; // defined in "filetypes/filetype_typedefs.h"
	int getLengthOfReal( ReinterpretReal reinterpret );
};

template <class FileType, class FilePattern, class MemoryPattern, class TheSite> LinkFile<FileType, FilePattern, MemoryPattern, TheSite>::LinkFile( ReinterpretReal reinterpret ) : filetype(  getLengthOfReal(reinterpret) ), reinterpret( reinterpret)
{
}

template <class FileType, class FilePattern, class MemoryPattern, class TheSite> LinkFile<FileType, FilePattern, MemoryPattern, TheSite>::~LinkFile()
{
}

/**
 *
 */
template <class FileType, class FilePattern, class MemoryPattern, class TheSite> bool LinkFile<FileType, FilePattern, MemoryPattern, TheSite>::load( TheSite site, std::string filename, Real *U )
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

	if( !filetype.loadHeader( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't read header");
		return false;
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

	if( reinterpret == STANDARD )
	{
		Real temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number
			file.read( (char*)&temp, sizeof(Real) );
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			if( index >= 0 ) U[index] = temp;
		}
	}
	else if( reinterpret == DOUBLE )
	{
		double temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number as DOUBLE
			file.read( (char*)&temp, sizeof(double) );
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			if( index >= 0 ) U[index] = (Real)temp; // cast to type Real
		}
	}
	else if( reinterpret == FLOAT )
	{
		float temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number AS FLOAT
			file.read( (char*)&temp, sizeof(float) );
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			if( index >= 0 ) U[index] = (Real)temp; // cast to type Real
		}
	}


	// load footer
	if( !filetype.loadFooter( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't read footer");
	}


	file.close();
	return true;
}


/**
 *
 */
template <class FileType, class FilePattern, class MemoryPattern, class TheSite> bool LinkFile<FileType, FilePattern, MemoryPattern, TheSite>::save( TheSite site, std::string filename, Real *U )
{
	// open file
	std::fstream file;
	file.open( filename.c_str(), std::ios::out | std::ios::binary);

	if( !file )
	{
		util::Logger::log( util::ERROR, "Can't open file");
		util::Logger::log( util::ERROR, filename.c_str() );
		return false;
	}

	// save header
	if( !filetype.saveHeader( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't write header");
		return false;
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







	if( reinterpret == STANDARD )
	{
		Real temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			temp = U[index];
			file.write( (char*)&temp, sizeof(Real) );
		}
	}
	else if( reinterpret == DOUBLE )
	{
		double temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			temp = (double)U[index];
			file.write( (char*)&temp, sizeof(double) );
		}
	}
	else if( reinterpret == FLOAT )
	{
		float temp;
		for( int i = 0; i < configSize; i++ )
		{
			// load number
			int index = MemoryPattern::getIndexByUnique( FilePattern::getUniqueIndex(i, site.size) , site.size );
			temp = (float)U[index];
			file.write( (char*)&temp, sizeof(float) );
		}
	}






	// save footer
	if( !filetype.saveFooter( &file ) )
	{
		util::Logger::log( util::ERROR, "Can't read footer");
	}


	file.close();
	return true;
}

// TODO place somewhere else
template <class FileType, class FilePattern, class MemoryPattern, class TheSite> int LinkFile<FileType, FilePattern, MemoryPattern, TheSite>::getLengthOfReal( ReinterpretReal reinterpret )
{
	if ( reinterpret == DOUBLE ) return sizeof( double );
	else if ( reinterpret == FLOAT ) return sizeof( float );
	else return sizeof( Real );
}


#endif /* LINKFILE_HXX_ */
