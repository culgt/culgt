/*
 * ProgramOptions.hxx
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
 */

#ifndef PROGRAMOPTIONS_HXX_
#define PROGRAMOPTIONS_HXX_

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>

#include <fstream>
#include <string>

using namespace std;

class ProgramOptions
{
public:
	ProgramOptions();
	int init( int argc, char* argv[] );

	string getFBasename() const {
		return fBasename;
	}

	string getFEnding() const {
		return fEnding;
	}

	int getFNumberformat() const {
		return fNumberformat;
	}

	int getFStartnumber() const {
		return fStartnumber;
	}

	FileType getFType() const {
		return fType;
	}

	int getGaugeCopies() const {
		return gaugeCopies;
	}

	int getNconf() const {
		return nconf;
	}

	bool isNoRandomTrafo() const {
		return noRandomTrafo;
	}

	int getOrCheckPrecision() const {
		return orCheckPrecision;
	}

	int getOrMaxIter() const {
		return orMaxIter;
	}

	float getOrParameter() const {
		return orParameter;
	}

	float getOrPrecision() const {
		return orPrecision;
	}

	ReinterpretReal getReinterpret() const {
		return reinterpret;
	}

	int getReproject() const {
		return reproject;
	}

	float getSaMax() const {
		return saMax;
	}

	int getSaMicroupdates() const {
		return saMicroupdates;
	}

	float getSaMin() const {
		return saMin;
	}

	int getSaSteps() const {
		return saSteps;
	}

	long getSeed() const {
		return seed;
	}

private:
	boost::program_options::variables_map options_vm;
	boost::program_options::options_description options_desc;

//	int argc;
//	char* argv[];

	// variables
	string configFile;

	FileType fType;
	string fBasename;
	string fEnding;
	int fNumberformat;
	int fStartnumber;
	int nconf;

	ReinterpretReal reinterpret;

	long seed;

	int gaugeCopies;
	bool noRandomTrafo;
	int reproject;

	int saSteps;
	float saMin;
	float saMax;
	int saMicroupdates;

	int orMaxIter;
	float orParameter;
	float orPrecision;
	int orCheckPrecision;



};

ProgramOptions::ProgramOptions() : options_desc("Allowed options")//, argc(argc), argv(*argv)
{
}

int ProgramOptions::init( int argc, char* argv[] )
{
	options_desc.add_options()
			("help", "produce help message")

			("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")

			("ftype", boost::program_options::value<FileType>(&fType), "type of configuration (PLAIN, HEADERONLY, VOGT)")
			("fbasename", boost::program_options::value<string>(&fBasename), "file basename (part before numbering starts)")
			("fending", boost::program_options::value<string>(&fEnding)->default_value(".vogt"), "file ending to append to basename (default: .vogt)")
			("fnumberformat", boost::program_options::value<int>(&fNumberformat)->default_value(1), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
			("fstartnumber", boost::program_options::value<int>(&fStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
			("nconf,m", boost::program_options::value<int>(&nconf)->default_value(1), "how many files to gaugefix")

			("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpret)->default_value(STANDARD), "reinterpret Real datatype (STANDARD = do nothing, FLOAT = read input as float and cast to Real, DOUBLE = ...)")

			("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")

			("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "Number of gauge copies")
			("norandomtrafo", boost::program_options::value<bool>(&noRandomTrafo)->default_value(false), "no random gauge trafo" )
			("reproject", boost::program_options::value<int>(&reproject)->default_value(100), "reproject every arg-th step")

			("sasteps", boost::program_options::value<int>(&saSteps)->default_value(1000), "number of SA steps")
			("samin", boost::program_options::value<float>(&saMin)->default_value(.01), "min. SA temperature")
			("samax", boost::program_options::value<float>(&saMax)->default_value(.4), "max. SA temperature")
			("microupdates", boost::program_options::value<int>(&saMicroupdates)->default_value(3), "number of microcanoncial updates at each SA temperature")

			("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(1000), "Max. number of OR iterations")
			("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "OR parameter")
			("orprecision", boost::program_options::value<float>(&orPrecision)->default_value(1E-7), "OR precision (dmuAmu)")
			("orcheckprecision", boost::program_options::value<int>(&orCheckPrecision)->default_value(100), "how often to check the gauge precision")
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
		cout << "Usage: " << argv[0] << " [options] [config-file]" << endl;
		cout << options_desc << "\n";
		return 1;
	}
	return 0;
}


#endif /* PROGRAMOPTIONS_HXX_ */
