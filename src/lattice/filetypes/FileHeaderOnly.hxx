/*
 * FileHeaderOnly.hxx
 * 
 * Class to read configurations in a "generic form" of the mdp-format,
 * i.e., the same data ordering as in the mdp-format is expected
 * but arbitrary headers are allowed (will simply be copied to the output file).
 * No footer allowed!
 *
 *  Created on: June 14, 2012
 *      Author: schroeck
 */

#ifndef FILEHEADERONLY_HXX_
#define FILEHEADERONLY_HXX_

#include <iostream>
#include <fstream>
//#include "../gaugefixing/GlobalConstants.hxx"

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
