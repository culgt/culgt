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
#include <iterator>

namespace culgt
{

class FileIterator
{
public:
	FileIterator(): fileBasename("out_"), fileExtension(".dat"), fileNumberformat(4), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileExtension ): fileBasename(fileBasename), fileExtension( fileExtension ), fileNumberformat(4), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileExtension, int fileNumberformat ): fileBasename(fileBasename), fileExtension( fileExtension ), fileNumberformat(fileNumberformat), fileNumberStart(0), fileNumberEnd(0), fileNumberStep(1), currentFileNumber(0)
	{
	};

	FileIterator( std::string fileBasename, std::string fileExtension, int fileNumberformat, int fileNumberStart, int fileNumberEnd, int fileNumberStep = 1 ): fileBasename(fileBasename), fileExtension( fileExtension ), fileNumberformat(fileNumberformat), fileNumberStart(fileNumberStart), fileNumberEnd(fileNumberEnd), fileNumberStep(fileNumberStep), currentFileNumber(fileNumberStart)
	{
	};

	std::string getFilename()
	{
		std::stringstream name;
		name << fileBasename << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << fileExtension;
		return name.str();
	}

	std::string getFilename( std::string appendix )
	{
		std::stringstream name;
		name << fileBasename << appendix << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << fileExtension;
		return name.str();
	}

	std::string getFilenameWithExtension( std::string extension  )
	{
		std::stringstream name;
		name << fileBasename << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << extension;
		return name.str();
	}

	std::string getFilenameWithExtension( std::string appendix, std::string extension  )
	{
		std::stringstream name;
		name << fileBasename << appendix << std::setw( fileNumberformat ) << std::setfill( '0' ) << currentFileNumber << extension;
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

	void setFileExtension(const std::string& fileExtension)
	{
		this->fileExtension = fileExtension;
	}

	class Iterator: public std::iterator<std::forward_iterator_tag,std::string>
	{
	private:
		int currentFileNumber;
		FileIterator& super;
	public:
		Iterator( int start, FileIterator& super ) :currentFileNumber( start ), super(super)
		{
		}
		std::string operator*()
		{
			super.currentFileNumber = currentFileNumber;
			return super.getFilename();
		}
		const Iterator* operator++()
		{
			currentFileNumber += super.fileNumberStep;
			return this;
		}
		bool operator==( const Iterator& rhs )
		{
			return currentFileNumber == rhs.currentFileNumber;
		}
		bool operator!=( const Iterator& rhs )
		{
			return currentFileNumber != rhs.currentFileNumber;
		}
	};

	Iterator begin()
	{
		return Iterator( Iterator( fileNumberStart, *this ) );
	}
	Iterator end()
	{
		return Iterator( Iterator( fileNumberEnd, *this ) );
	}

private:
	std::string fileBasename;
	std::string fileExtension;

	int fileNumberformat;

	int fileNumberStart;
	int fileNumberEnd;
	int fileNumberStep;

	int currentFileNumber;
};

}



#endif /* FILEITERATOR_H_ */
