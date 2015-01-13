/*
 * AutotuneManager.h
 *
 *  Created on: Oct 30, 2014
 *      Author: vogt
 */

#ifndef AUTOTUNEMANAGER_H_
#define AUTOTUNEMANAGER_H_

#include <string>
#include <sstream>
#include <fstream>
#include <boost/functional/hash.hpp>

namespace culgt
{

class AutotuneManagerOptionNotAvailable: public std::exception
{
public:
	~AutotuneManagerOptionNotAvailable() throw() {};
	AutotuneManagerOptionNotAvailable(){};
};


class AutotuneManager
{
public:
	AutotuneManager()
	{
		identifier << "tune";
	}

	void addHashedAttribute( std::string attribute )
	{
		boost::hash<std::string> string_hash;
		std::size_t h = string_hash(attribute);
		identifier << "_"  << h;
	}

	void addAttribute( std::string attribute )
	{
		identifier << "_" << attribute;
	}

	void addAttribute( int attribute )
	{
		identifier << "_" << attribute;
	}

	std::string getIdentifier()
	{
		return identifier.str();
	}

	RuntimeChooserOption getOptimalId()
	{
		std::ifstream file( getIdentifier().c_str() );
		RuntimeChooserOption optimalId;
		if( file.good() )
		{
			file >> optimalId.name;
			file >> optimalId.id;
		}
		else
		{
			throw AutotuneManagerOptionNotAvailable();
		}
		file.close();
		return optimalId;
	}

	void writeOptimalId( RuntimeChooserOption optimalId )
	{
		std::ofstream file( getIdentifier().c_str() );
		file << optimalId.name << endl;
		file << optimalId.id;
		file.close();
	}

private:
	std::stringstream identifier;
};

}
#endif /* AUTOTUNEMANAGER_H_ */
