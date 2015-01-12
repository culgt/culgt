/**
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef TUNABLEOBJECT_H_
#define TUNABLEOBJECT_H_

#include "../../cuLGT1legacy/Chronotimer.h"
#include "AutotuneManager.h"

#include "cuda.h"
#include <vector>

using std::vector;

namespace culgt
{

class LastElementReachedException: public std::exception
{
public:
	~LastElementReachedException() throw() {};
	LastElementReachedException(){};
};

class InvalidKernelSetupException: public std::exception
{
public:
	~InvalidKernelSetupException() throw() {};
	InvalidKernelSetupException(){};
};


class TunableObject
{
public:
	TunableObject(){};
	virtual ~TunableObject(){};
	virtual void run( size_t id ) = 0;
	virtual void run() = 0;
	virtual void preTuneAction() = 0;
	virtual std::vector<size_t>* getOptions() = 0;

	virtual string getClassName()
	{
		return typeid( *this ).name();
	}

	double measure( size_t id, int iter = 100 )
	{
		timer.reset();
		timer.start();
		CUDA_LAST_ERROR("TunableObject::reset()/start()" );
		cudaDeviceSynchronize();

		for( int i = 0; i < iter; i++ )
		{
			run( id );
		}

		cudaDeviceSynchronize();
		timer.stop();
		CUDA_LAST_ERROR("TunableObject::stop()" );
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

	void tune( int iter = 1000, bool force = false )
	{
		AutotuneManager autotuneManager;
		autotuneManager.addHashedAttribute( getClassName() );
		int selectedDeviceNumber;
		cudaGetDevice( &selectedDeviceNumber );
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop,selectedDeviceNumber);
		autotuneManager.addAttribute( string(prop.name) );

		try
		{
			optimalId = autotuneManager.getOptimalId();
		}
		catch( AutotuneManagerOptionNotAvailable& e )
		{
			force = true;
		}

		if( force )
		{
			forceTune( iter );
			autotuneManager.writeOptimalId( optimalId );
		}
		std::cout << "Using option: " << optimalId << std::endl;
	}

	void forceTune( int iter )
	{
		double bestPerformance = 0;

		std::vector<size_t>* ptrToOptions = getOptions();

		std::cout << ptrToOptions->size() << std::endl;

		for( vector<size_t>::iterator it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
		{
			try
			{
				std::cout << "warmup option " << *it << std::endl;
				measure( *it, iter );
			}
			catch(InvalidKernelSetupException& e)
			{
				continue;
			}
			std::cout << "warmup done..." << std::endl;
			break;
		}

		for( vector<size_t>::iterator it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
		{
			preTuneAction();

			double performance;
			try
			{
				performance = 1./measure( *it, iter );
			}
			catch(InvalidKernelSetupException& e)
			{
				std::cout << "Tuning option " << *it << ": invalid configuration" << std::endl;
//				id++;
				continue;
			}
			CUDA_LAST_ERROR("TunableObject::measure()" );

			std::cout << "Tuning option " << *it << ": " << performance << std::endl;


			if( performance > bestPerformance )
			{
				bestPerformance = performance;
				optimalId = *it;
			}
		}
	}

protected:
	size_t optimalId;

private:
	Chronotimer timer;
};

}

#endif /* MEASURABLEOBJECT_H_ */
