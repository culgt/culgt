/*
 * HeaderVogt.hxx
 *
 *  Created on: Apr 30, 2012
 *      Author: vogt
 */

#ifndef HEADERVOGT_HXX_
#define HEADERVOGT_HXX_

#include <fstream>

class HeaderVogt
{
public:
	HeaderVogt();
	virtual ~HeaderVogt();
	bool load( std::fstream* file );
};


HeaderVogt::HeaderVogt()
{
}

HeaderVogt::~HeaderVogt()
{
}


bool HeaderVogt::load( std::fstream* file )
{
	short tempShort;


	std::cout << "size of short: " << sizeof(short) <<  std::endl;

	file->read( (char*)&tempShort, sizeof(short) ); // ndim
	std::cout << tempShort << std::endl;
	file->read( (char*)&tempShort, sizeof(short) ); // nc
	std::cout << tempShort << std::endl;

	file->read( (char*)&tempShort, sizeof(short) ); // t
	std::cout << tempShort << std::endl;
	file->read( (char*)&tempShort, sizeof(short) ); // x
	std::cout << tempShort << std::endl;
	file->read( (char*)&tempShort, sizeof(short) ); // y
	std::cout << tempShort << std::endl;
	file->read( (char*)&tempShort, sizeof(short) ); // z
	std::cout << tempShort << std::endl;

	file->read( (char*)&tempShort, sizeof(short) ); // length of real
	std::cout << tempShort << std::endl;





	return true;
}

#endif /* HEADERVOGT_HXX_ */
