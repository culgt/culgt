/**
 *
 *  Created on: Apr 14, 2014
 *      Author: vogt
 */

#ifndef FILEITERATOR_H_
#define FILEITERATOR_H_


#include <string>
#include <sstream>
#include <iomanip>

namespace culgt
{

class FileIterator
{
public:
	FileIterator(): fileBasename("out_"), fileEnding(".dat"), fileNumberformat(4), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileEnding ): fileBasename(fileBasename), fileEnding( fileEnding ), fileNumberformat(4), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileEnding, int fileNumberformat ): fileBasename(fileBasename), fileEnding( fileEnding ), fileNumberformat(fileNumberformat), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileEnding, int fileNumberformat, int fileNumberStart, int fileNumberEnd, int fileNumberStep = 1 ): fileBasename(fileBasename), fileEnding( fileEnding ), fileNumberformat(fileNumberformat), fileNumberStart(fileNumberStart), fileNumberEnd(fileNumberEnd), fileNumberStep(fileNumberStep), currentFileNumber(fileNumberStart)
	{
	};

	std::string getFilename()
	{
		std::stringstream name;
		name << fileBasename << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << fileEnding;
		return name.str();
	}

	std::string getFilename( std::string appendix )
	{
		std::stringstream name;
		name << fileBasename << appendix << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << fileEnding;
		return name.str();
	}

	bool hasElement()
	{
		return currentFileNumber <= fileNumberEnd;
	}

	std::string next()
	{
		std::string result = getFilename();
		currentFileNumber += fileNumberStep;
		return result;
	}

	void reset()
	{
		currentFileNumber = fileNumberStart;
	}
private:
	std::string fileBasename;
	std::string fileEnding;

	int fileNumberformat;

	int fileNumberStart;
	int fileNumberEnd;
	int fileNumberStep;

	int currentFileNumber;
};

}



#endif /* FILEITERATOR_H_ */
