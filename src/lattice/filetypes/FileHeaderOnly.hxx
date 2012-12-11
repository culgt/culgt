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
 * Class to read configurations in a "generic form" of the mdp-format,
 * i.e., the same data ordering as in the mdp-format is expected
 * but arbitrary headers are allowed (will simply be copied to the output file).
 * No footer allowed!
 *
 * TODO give an explict definition of th format here!
 *
 *  Created on: June 14, 2012
 *      Author: schroeck
 */

#ifndef FILEHEADERONLY_HXX_
#define FILEHEADERONLY_HXX_

#include <iostream>
#include <fstream>
#include "../../gaugefixing/GlobalConstants.h" // TODO must not use a gaugefixing class!!!

using namespace std;

class FileHeaderOnly
{
public:
	FileHeaderOnly( int LENGTH_OF_REAL );
	virtual ~FileHeaderOnly();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	bool loadFooter( std::fstream* file );
	bool saveFooter( std::fstream* file );
private:
	long arraySize;
	long filelength;
	int offset;
	char *header;
	int LENGTH_OF_REAL;
};


FileHeaderOnly::FileHeaderOnly( int LENGTH_OF_REAL ) : LENGTH_OF_REAL( LENGTH_OF_REAL )
{
	int Ndim = 4;
	int Nc = 3;
	arraySize=Ndim*Nc*Nc*2*Nx*Ny*Nz*Nt*LENGTH_OF_REAL;
}

FileHeaderOnly::~FileHeaderOnly()
{
}


bool FileHeaderOnly::loadHeader( std::fstream* file )
{
	//what's the offset (the header length)?
	file->seekg(0, ios::end);
  filelength = file->tellg();
  file->seekg(0, ios::beg);
	offset = filelength - arraySize;
	cout << "offset by header = " << offset << '\n' << endl;
	
	header = (char*) malloc(offset);

  //copy header
  file->read(header,offset);
	
	if( file->fail() )
		return false;

	cout << header << endl;

	return true;
}


bool FileHeaderOnly::saveHeader( std::fstream* file )
{
	file->write(header,offset);
	
	if( file->fail() )
		return false;
	else
		return true;
}

bool FileHeaderOnly::loadFooter( std::fstream* file )
{
	return true;
}

bool FileHeaderOnly::saveFooter( std::fstream* file )
{
	return true;
}


#endif /* FILEMDP_HXX_ */
