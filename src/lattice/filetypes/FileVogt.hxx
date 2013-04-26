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
 * TODO give definition of filetype here!
 *
 * TODO: MAKE A LATTICE SIZE CHECK!!!!!
 *
 */

#ifndef FILEVOGT_HXX_
#define FILEVOGT_HXX_

#include <fstream>

class FileVogt
{
public:
	FileVogt( int LENGTH_OF_REAL );
	virtual ~FileVogt();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	bool loadFooter( std::fstream* file );
	bool saveFooter( std::fstream* file );
	short ndim;
	short nc;
	short latsize[20];
	short lengthOfReal;
private:
	int LENGTH_OF_REAL;
};


FileVogt::FileVogt( int LENGTH_OF_REAL ) : LENGTH_OF_REAL( LENGTH_OF_REAL )
{
}

FileVogt::~FileVogt()
{
}


bool FileVogt::loadHeader( std::fstream* file )
{
//	short tempShort;


	std::cout << "size of short: " << sizeof(short) <<  std::endl;

	file->read( (char*)&ndim, sizeof(short) ); // ndim
	std::cout << ndim << std::endl;
	file->read( (char*)&nc, sizeof(short) ); // nc
	std::cout << nc << std::endl;

	for( int i = 0; i < ndim; i++ )
	{
		file->read( (char*)&latsize[i], sizeof(short) );
		std::cout << latsize[i] << std::endl;
	}

	file->read( (char*)&lengthOfReal, sizeof(short) ); // length of real
	std::cout << lengthOfReal << std::endl;

	if( lengthOfReal != LENGTH_OF_REAL ) std::cout << "WRONG LENGTH OF REAL: " << lengthOfReal << " bytes in header, while app wants " << LENGTH_OF_REAL << " bytes" << std::endl;

	return true;
}

bool FileVogt::saveHeader( std::fstream* file )
{
//	short tempShort;


	std::cout << "size of short: " << sizeof(short) <<  std::endl;

	file->write( (char*)&ndim, sizeof(short) ); // ndim
	std::cout << ndim << std::endl;
	file->write( (char*)&nc, sizeof(short) ); // nc
	std::cout << nc << std::endl;

	for( int i = 0; i < ndim; i++ )
	{
		file->write( (char*)&latsize[i], sizeof(short) );
		std::cout << latsize[i] << std::endl;
	}

	file->write( (char*)&lengthOfReal, sizeof(short) ); // length of real
	std::cout << lengthOfReal << std::endl;

	return true;
}

bool FileVogt::loadFooter( std::fstream* file )
{
	return true;
}

bool FileVogt::saveFooter( std::fstream* file )
{
	return true;
}


#endif /* FILEVOGT_HXX_ */
