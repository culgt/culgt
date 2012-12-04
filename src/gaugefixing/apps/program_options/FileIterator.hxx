/*
 * FileIterator.hxx
 *
 * Use a std-lib like iterator...
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
