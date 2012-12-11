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
 * Same format as HeaderOnly, but no header allowed => obsolete...
 *
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
