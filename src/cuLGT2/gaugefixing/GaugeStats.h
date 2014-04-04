/**
 *
 *  Created on: Mar 26, 2014
 *      Author: vogt
 */

#ifndef GAUGESTATS_H_
#define GAUGESTATS_H_

class GaugeStats
{
public:
	GaugeStats(): gff(0), precision(0)
	{
	}
	GaugeStats( double gff, double precision ): gff(gff), precision(precision)
	{
	}

	double getGff() const
	{
		return gff;
	}

	double getPrecision() const
	{
		return precision;
	}

private:
	double gff;
	double precision;
};


#endif /* GAUGESTATS_H_ */
