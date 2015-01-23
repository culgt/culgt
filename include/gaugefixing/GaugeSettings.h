/**
 * GaugeSettings.h
 *
 *  Created on: Apr 23, 2014
 *      Author: vogt
 */

#ifndef GAUGESETTINGS_H_
#define GAUGESETTINGS_H_

#include <boost/program_options/options_description.hpp>

namespace culgt
{

class GaugeSettings
{
public:
	int getCheckPrecision() const
	{
		return checkPrecision;
	}

	void setCheckPrecision(int checkPrecision)
	{
		this->checkPrecision = checkPrecision;
	}

	int getGaugeCopies() const
	{
		return gaugeCopies;
	}

	void setGaugeCopies(int gaugeCopies)
	{
		this->gaugeCopies = gaugeCopies;
	}

	bool isRandomTrafo() const
	{
		return randomTrafo;
	}

	void setRandomTrafo(bool randomTrafo)
	{
		this->randomTrafo = randomTrafo;
	}

	int getOrMaxIter() const
	{
		return orMaxIter;
	}

	void setOrMaxIter(int orMaxIter)
	{
		this->orMaxIter = orMaxIter;
	}

	float getOrParameter() const
	{
		return orParameter;
	}

	void setOrParameter(float orParameter)
	{
		this->orParameter = orParameter;
	}

	double getPrecision() const
	{
		return precision;
	}

	void setPrecision(double precision)
	{
		this->precision = precision;
	}

	bool isPrintStats() const
	{
		return printStats;
	}

	void setPrintStats(bool printStats)
	{
		this->printStats = printStats;
	}

	int getReproject() const
	{
		return reproject;
	}

	void setReproject(int reproject)
	{
		this->reproject = reproject;
	}

	float getSaMax() const
	{
		return saMax;
	}

	void setSaMax(float saMax)
	{
		this->saMax = saMax;
	}

	float getSaMin() const
	{
		return saMin;
	}

	void setSaMin(float saMin)
	{
		this->saMin = saMin;
	}

	int getSaSteps() const
	{
		return saSteps;
	}

	void setSaSteps(int saSteps)
	{
		this->saSteps = saSteps;
	}

	int getMicroiter() const
	{
		return microiter;
	}

	void setMicroiter( int microiter )
	{
		this->microiter = microiter;
	}

//	long getSeed() const
//	{
//		return seed;
//	}
//
//	void setSeed(long seed)
//	{
//		this->seed = seed;
//	}

//	void setGaugeOptions( boost::program_options::options_description* po )
	boost::program_options::options_description getGaugeOptions()
	{
		boost::program_options::options_description gaugeOptions("Gaugefixing options");

//		po->add_options()
		gaugeOptions.add_options()
//					("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")

					("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "calculate <arg> gauge copies and save the one with best functional value")
					("randomtrafo", boost::program_options::value<bool>(&randomTrafo)->default_value(true), "apply a random gauge transformation before gauge fixing")

					("precision", boost::program_options::value<double>(&precision)->default_value(1e-14), "desired precision")
					("checkprecision", boost::program_options::value<int>(&checkPrecision)->default_value(100), "check if precision is reached every <arg> step")
					("reproject", boost::program_options::value<int>(&reproject)->default_value(100), "reproject links to SU(N) every <arg> step")

					("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(100), "maximal number of overrelaxation iterations")
					("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "overrelaxation parameter w (g^w)")

					("samax", boost::program_options::value<float>(&saMax)->default_value(1.4), "Simulated Annealing start temperature")
					("samin", boost::program_options::value<float>(&saMin)->default_value(0.1), "Simulated Annealing end temperature")
					("sasteps", boost::program_options::value<int>(&saSteps)->default_value(0), "number of Simulated Annealing steps")

					("microiter", boost::program_options::value<int>(&microiter)->default_value(3), "number of microcanonical updates per heatbath in Simulated Annealing")

					("printstats", boost::program_options::value<bool>(&printStats)->default_value(true), "print progress on command line");

		return gaugeOptions;
	}


private:

	double precision;

	int saSteps;
	float saMax;
	float saMin;

	int microiter;

	int orMaxIter;
	float orParameter;

	int checkPrecision;
	int reproject;

	int gaugeCopies;

//	long seed;

	bool randomTrafo;

	bool printStats;
};

}

#endif /* GAUGESETTINGS_H_ */
