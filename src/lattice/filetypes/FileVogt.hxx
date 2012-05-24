/*
 * HeaderVogt.hxx
 *
 *  Created on: Apr 30, 2012
 *      Author: vogt
 */

#ifndef FILEVOGT_HXX_
#define FILEVOGT_HXX_

#include <fstream>

class FileVogt
{
public:
	FileVogt();
	virtual ~FileVogt();
	bool loadHeader( std::fstream* file );
	bool saveHeader( std::fstream* file );
	short ndim;
	short nc;
	short latsize[20];
	short lengthOfReal;
};


FileVogt::FileVogt()
{
}

FileVogt::~FileVogt()
{
}


bool FileVogt::loadHeader( std::fstream* file )
{
//	short tempShort;


	std::cout << "size of short: " << sizeof(short) <<  std::endl;

	file->read( (char*)&ndim, sizeof(short) ); // ndim
	std::cout << ndim << std::endl;
	file->read( (char*)&nc, sizeof(short) ); // nc
	std::cout << nc << std::endl;

	for( int i = 0; i < ndim; i++ )
	{
		file->read( (char*)&latsize[i], sizeof(short) );
		std::cout << latsize[i] << std::endl;
	}

	file->read( (char*)&lengthOfReal, sizeof(short) ); // length of real
	std::cout << lengthOfReal << std::endl;

	return true;
}

bool FileVogt::saveHeader( std::fstream* file )
{
//	short tempShort;


	std::cout << "size of short: " << sizeof(short) <<  std::endl;

	file->write( (char*)&ndim, sizeof(short) ); // ndim
	std::cout << ndim << std::endl;
	file->write( (char*)&nc, sizeof(short) ); // nc
	std::cout << nc << std::endl;

	for( int i = 0; i < ndim; i++ )
	{
		file->write( (char*)&latsize[i], sizeof(short) );
		std::cout << latsize[i] << std::endl;
	}

	file->write( (char*)&lengthOfReal, sizeof(short) ); // length of real
	std::cout << lengthOfReal << std::endl;

	return true;
}


#endif /* FILEVOGT_HXX_ */
