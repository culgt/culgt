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
#include "../template_instantiation/RuntimeChooserOption.h"

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
		attribute.erase(remove_if(attribute.begin(), attribute.end(), ::isspace), attribute.end());
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
		if( optimalId.name.size() == 0 )
			file << "N/A\n";
		else
			file << optimalId.name << "\n";
		file << optimalId.id << std::endl;
		file.close();
	}

private:
	std::stringstream identifier;
};

}
#endif /* AUTOTUNEMANAGER_H_ */
