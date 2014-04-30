/**
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef TUNABLEOBJECT_H_
#define TUNABLEOBJECT_H_

#include "../../cuLGT1legacy/Chronotimer.h"

namespace culgt
{

class LastElementReachedException: public std::exception
{
public:
	~LastElementReachedException() throw() {};
	LastElementReachedException(){};
};


class TunableObject
{
public:
	TunableObject(){};
	virtual ~TunableObject(){};
	virtual void run( int id ) = 0;
	virtual void preTuneAction() = 0;
	double measure( int id, int iter = 100 )
	{
		timer.reset();
		timer.start();
		cudaDeviceSynchronize();

		for( int i = 0; i < iter; i++ )
		{
			run( id );
		}

		cudaDeviceSynchronize();
		timer.stop();
		return timer.getTime();
	}

	void startTime()
	{
		timer.reset();
		timer.start();
		cudaDeviceSynchronize();
	}

	void stopTime()
	{
		cudaDeviceSynchronize();
		timer.stop();
	}

	double getTime()
	{
		cudaDeviceSynchronize();
		return timer.getTime();
	}

	void tune( int iter = 1000 )
	{
		double bestPerformance = 0;
		int id = 0;

		// warmup
		measure( 0 );

		while( true )
		{
			try
			{
				preTuneAction();
				double performance = 1./measure( id, iter );

				std::cout << "Tuning option " << id << ": " << performance << std::endl;


				if( performance > bestPerformance )
				{
					bestPerformance = performance;
					optimalId = id;
				}
			}
			catch( LastElementReachedException& e )
			{
				break;
			}
			id++;
		}
		std::cout << "Using option: " << optimalId << std::endl;
	}

protected:
	int optimalId;

private:
	Chronotimer timer;
};

}

#endif /* MEASURABLEOBJECT_H_ */
