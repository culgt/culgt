/**
 * GaugeSettings.h
 *
 *  Created on: Apr 23, 2014
 *      Author: vogt
 */

#ifndef GAUGESETTINGS_H_
#define GAUGESETTINGS_H_

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

	int getSeed() const
	{
		return seed;
	}

	void setSeed(int seed)
	{
		this->seed = seed;
	}

private:
	double precision;

	int saSteps;
	float saMax;
	float saMin;

	int orMaxIter;
	float orParameter;

	int checkPrecision;
	int reproject;

	int gaugeCopies;

	int seed;

	bool randomTrafo;

	bool printStats;
};

}

#endif /* GAUGESETTINGS_H_ */
