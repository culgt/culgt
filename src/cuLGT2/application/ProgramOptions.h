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

#include "../lattice/filetype_typedefs.h"

using std::string;
using std::ifstream;

namespace culgt
{



class ProgramOptions
{
public:
	ProgramOptions()
	{
		options_desc.add_options()
		("help", "produce help message")

		("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")

		("nt", boost::program_options::value<int>(&nt)->default_value(32), "size in t direction" )
		("nx", boost::program_options::value<int>(&nx)->default_value(32), "size in x direction" )
		("ny", boost::program_options::value<int>(&ny)->default_value(-1), "size in y direction (-1 = same as x direction)" )
		("nz", boost::program_options::value<int>(&nz)->default_value(-1), "size in z direction (-1 = same as x direction)" )

		("fbasename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
		("fending", boost::program_options::value<string>(&fileEnding)->default_value(".dat"), "file ending to append to basename (default: .vogt)")
		("fnumberformat", boost::program_options::value<int>(&fileNumberformat)->default_value(4), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
		("fstartnumber", boost::program_options::value<int>(&fileNumberStart)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
		("fstepnumber", boost::program_options::value<int>(&fileNumberStep)->default_value(1), "load every <fstepnumber>-th file")
		("nconf,m", boost::program_options::value<int>(&nConf)->default_value(1), "how many files to iterate")
		("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpretReal)->default_value(STANDARD), "interprets the real values in file as <arg>")


		("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")
		;
	}

#if BOOST_VERSION < 105300
	void parseOptions( int argc, char* argv[], bool beQuiet = false )
#else
	void parseOptions( const int argc, const char* const argv[], bool beQuiet = false )
#endif
	{
		try
			{

			if( beQuiet )
			{

	//			boost::program_options::command_line_parser clp(argc, argv);
	////			clp.options(options_desc).allow_unregistered();
	//			clp.positional(options_p).options(options_desc.add(hiddenCatchAll)).allow_unregistered();
	//
	//			boost::program_options::parsed_options parsed = clp.run();
	//
	//			boost::program_options::store( parsed, options_vm );
	//
	//			std::vector<std::string> unrecognized = boost::program_options::collect_unrecognized( parsed.options, boost::program_options::exclude_positional );
	//
	//			if (!unrecognized.empty())
	//			{
	//				std::cout << "Got " << unrecognized.size() << " unrecognized option(s):" << std::endl;
	//				for( int i = 0; i < unrecognized.size(); i++ )
	//					std::cout << '\t' << unrecognized[i] << std::endl;
	//			}
	//
	////			if (!catchAll.empty())
	//			{
	//				std::cout << "Got " << catchAll.size() << " catchAll option(s):" << std::endl;
	//				for( int i = 0; i < catchAll.size(); i++ )
	//					std::cout << '\t' << catchAll[i] << std::endl;
	//			}
	//


				boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
						options(options_desc).allow_unregistered().run(), options_vm);
			}
			else
			{
				boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
						options(options_desc).run(), options_vm);
			}
			boost::program_options::notify(options_vm);

			if( configFile.length() > 0 )
			{
				ifstream cfg( configFile.c_str() );
				if( !cfg )
				{
					showError( "Config file does not exist!" );
				}

				boost::program_options::store(boost::program_options::parse_config_file( cfg, options_desc, beQuiet), options_vm);
				boost::program_options::notify(options_vm);
			}

			if( !beQuiet )
			{
				if (options_vm.count("help")) {

					showOptionsAndExit( argv[0] );
				}
			}
		}
		catch( boost::program_options::error& e )
		{
			std::cout <<  "\n" << e.what() << "\n"<< std::endl;

			showOptionsAndExit( argv[0] );
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

	void addOption( boost::program_options::options_description option )
	{
		options_desc.add( option );
	}

//	boost::program_options::options_description* getOptionsDescription()
//	{
//		return &options_desc;
//	}

	long getSeed() const
	{
		return seed;
	}

	ReinterpretReal getReinterpretReal() const
	{
		return reinterpretReal;
	}

	int getNt() const
	{
		return nt;
	}

	int getNx() const
	{
		return nx;
	}

	int getNy() const
	{
		if( ny == -1 )
			return nx;
		else
			return ny;
	}

	int getNz() const
	{
		if( nz == -1 )
			return nx;
		else
			return nz;
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

	int nt;
	int nx;
	int ny;
	int nz;

	long seed;

	ReinterpretReal reinterpretReal;

	void showError( std::string message )
	{
		std::cout << message << std::endl;
		exit(1);
	}

	void showOptionsAndExit( const char* programName )
	{
		std::cout << "Usage: " << programName << " [options]" << std::endl;
		std::cout << options_desc << "\n";
		exit(1);
	}
};




}
#endif /* PROGRAMOPTIONS_H_ */
