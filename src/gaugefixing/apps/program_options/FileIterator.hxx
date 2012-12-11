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
 *  TODO Use a std-lib like iterator...
 */

#ifndef FILEITERATOR_HXX_
#define FILEITERATOR_HXX_

#include "ProgramOptions.hxx"
#include <iomanip>
#include <sstream>


class FileIterator
{
public:
	FileIterator( ProgramOptions options );
	bool hasNext();
	bool next();
	void reset();
	string getFilename();
	string getOutputFilename();
private:
	int currentId;
	ProgramOptions options;
};


FileIterator::FileIterator( ProgramOptions options ) : options(options)
{
	reset();
}

void FileIterator::reset()
{
	currentId = options.getFStartnumber();
}

bool FileIterator::hasNext()
{
	if( currentId < options.getFStartnumber()+options.getNconf()*options.getFStepnumber() )
		return true;
	else
		return false;
}

bool FileIterator::next()
{
	currentId += options.getFStepnumber();
	if( currentId  < options.getFStartnumber()+options.getNconf()*options.getFStepnumber() )
		return true;
	else
		return false;
}

string FileIterator::getFilename()
{
	stringstream filename(stringstream::out);
	filename << options.getFBasename() << setw( options.getFNumberformat() ) << setfill( '0' ) << currentId << options.getFEnding();
	return filename.str();
}

string FileIterator::getOutputFilename()
{
	stringstream filename(stringstream::out);
	filename << options.getFBasename() << options.getFOutputAppendix() << setw( options.getFNumberformat() ) << setfill( '0' ) << currentId << options.getFEnding();
	return filename.str();
}


#endif /* FILEITERATOR_HXX_ */
