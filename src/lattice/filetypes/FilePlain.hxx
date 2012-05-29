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
	FilePlain();
	virtual ~FilePlain();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	bool loadFooter( std::fstream* file );
	bool saveFooter( std::fstream* file );
	short ndim;
	short nc;
	short latsize[20];
	short lengthOfReal;
};


FilePlain::FilePlain()
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
