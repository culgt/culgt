/**
 * ProgramOptions.h
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef PROGRAMOPTIONS_H_
#define PROGRAMOPTIONS_H_

#include <string>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <iosfwd>

using std::string;
using std::ifstream;

namespace culgt
{



class ProgramOptions
{
public:
	ProgramOptions( const int argc, const char* argv[] )
	{
		options_desc.add_options()
		("help", "produce help message")

		("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")

		("fbasename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
		("fending", boost::program_options::value<string>(&fileEnding)->default_value(".dat"), "file ending to append to basename (default: .vogt)")
		("fnumberformat", boost::program_options::value<int>(&fileNumberformat)->default_value(4), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
		("fstartnumber", boost::program_options::value<int>(&fileNumberStart)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
		("fstepnumber", boost::program_options::value<int>(&fileNumberStep)->default_value(1), "load every <fstepnumber>-th file")
		("nconf,m", boost::program_options::value<int>(&nConf)->default_value(1), "how many files to iterate")
		;

		boost::program_options::positional_options_description options_p;
		options_p.add("config-file", -1);

		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(options_desc).positional(options_p).run(), options_vm);
		boost::program_options::notify(options_vm);

		ifstream cfg( configFile.c_str() );
		boost::program_options::store(boost::program_options::parse_config_file( cfg, options_desc), options_vm);
		boost::program_options::notify(options_vm);

		if (options_vm.count("help")) {
			std::cout << "Usage: " << argv[0] << " [options] [config-file]" << std::endl;
			std::cout << options_desc << "\n";
			exit(1);
		}
	}

	const string& getConfigFile() const
	{
		return configFile;
	}

	const string& getFileBasename() const
	{
		return fileBasename;
	}

	const string& getFileEnding() const
	{
		return fileEnding;
	}

	int getFileNumberformat() const
	{
		return fileNumberformat;
	}

	int getFileNumberStart() const
	{
		return fileNumberStart;
	}

	int getFileNumberStep() const
	{
		return fileNumberStep;
	}

	int getNConf() const
	{
		return nConf;
	}

private:
	boost::program_options::variables_map options_vm;
	boost::program_options::options_description options_desc;

	string configFile;

	string fileBasename;
	string fileEnding;

	int fileNumberformat;

	int fileNumberStart;
	int fileNumberStep;
	int nConf;
};




}
#endif /* PROGRAMOPTIONS_H_ */
