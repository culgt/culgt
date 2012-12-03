/*
 * HeaderVogt.hxx
 *
 *  Created on: Apr 30, 2012
 *      Author: vogt
 */

#ifndef FILEPLAIN_HXX_
#define FILEPLAIN_HXX_

#include <fstream>

class FilePlain
{
public:
	FilePlain( int LENGTH_OF_REAL );
	virtual ~FilePlain();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	bool loadFooter( std::fstream* file );
	bool saveFooter( std::fstream* file );
private:
	int LENGTH_OF_REAL;
};


FilePlain::FilePlain( int LENGTH_OF_REAL ) : LENGTH_OF_REAL(LENGTH_OF_REAL)
{
}

FilePlain::~FilePlain()
{
}


bool FilePlain::loadHeader( std::fstream* file )
{
	return true;
}

bool FilePlain::saveHeader( std::fstream* file )
{
	return true;
}

bool FilePlain::loadFooter( std::fstream* file )
{
	return true;
}

bool FilePlain::saveFooter( std::fstream* file )
{
	return true;
}


#endif /* FILEPLAIN_HXX_ */
