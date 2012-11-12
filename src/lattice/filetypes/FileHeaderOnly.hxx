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
#include "../gaugefixing/GlobalConstants.hxx"

using namespace std;

class FileHeaderOnly
{
public:
	FileHeaderOnly();
	virtual ~FileHeaderOnly();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	bool loadFooter( std::fstream* file );
	bool saveFooter( std::fstream* file );
	short ndim;
	short nc;
	short latsize[20];
	short lengthOfReal;
	long arraySize;
private:
	long filelength;
	int offset;
	char *header;
};


FileHeaderOnly::FileHeaderOnly()
{
//	//TODO get the lattice constants from main
	const lat_dim_t Ndim = 4;
	const short Nc = 3;
//	#ifdef _X_
//	const lat_coord_t Nx = _X_;
//	#else
//	#error "Define X (the lattice size in x-direction)"
//	#endif
//	#ifdef _Y_
//	const lat_coord_t Ny = _Y_;
//	#else
//	const lat_coord_t Ny = _X_;
//	bool warnY = true;
//	#endif
//	#ifdef _Z_
//	const lat_coord_t Nz = _Z_;
//	#else
//	const lat_coord_t Nz = _X_;
//	bool warnZ = true;
//	#endif
//	#ifdef _T_
//	const lat_coord_t Nt = _T_;
//	#else
//	#error "Define T (the lattice size in t-direction)"
//	#endif
//	//-----------------
	
#ifdef FILESP
	arraySize=Ndim*Nc*Nc*2*Nx*Ny*Nz*Nt*sizeof(float);
#else
	arraySize=Ndim*Nc*Nc*2*Nx*Ny*Nz*Nt*sizeof(double);
#endif

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
