/**
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef TUNABLEOBJECT_H_
#define TUNABLEOBJECT_H_

#include <vector>
#include "util/timer/Chronotimer.h"
#include "../template_instantiation/SequenceRunner.h"
#include "../template_instantiation/RuntimeChooser.h"
#include "AutotuneManager.h"

#include "cudacommon/cuda_error.h"
#include "cudacommon/DeviceProperties.h"

namespace culgt
{

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
	virtual std::vector<RuntimeChooserOption>* getOptions() = 0;

	virtual string getClassName()
	{
		return typeid( *this ).name();
	}

	double measure( size_t id, int iter = 100 )
	{
		timer.reset();
		timer.start();
		CUDA_LAST_ERROR("TunableObject::reset()/start()" );
#ifdef __CUDACC__
		cudaDeviceSynchronize();
#endif

		for( int i = 0; i < iter; i++ )
		{
			run( id );
		}

#ifdef __CUDACC__
		cudaDeviceSynchronize();
#endif
		timer.stop();
		CUDA_LAST_ERROR("TunableObject::stop()" );
		return timer.getTime();
	}

	void startTime()
	{
		timer.reset();
		timer.start();
#ifdef __CUDACC__
		cudaDeviceSynchronize();
#endif
	}

	void stopTime()
	{
#ifdef __CUDACC__
		cudaDeviceSynchronize();
#endif
		timer.stop();
	}

	double getTime()
	{
#ifdef __CUDACC__
		cudaDeviceSynchronize();
#endif
		return timer.getTime();
	}

	void tune( int iter = 1000, bool force = false )
	{
		AutotuneManager autotuneManager;
		autotuneManager.addHashedAttribute( getClassName() );
		autotuneManager.addAttribute( DeviceProperties::getName() );

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
			if( forceTune( iter ) )
				autotuneManager.writeOptimalId( optimalId );
			else
			{
				std::cout << "Autotune failed! Exiting..." << std::endl;
				exit(-1);
			}
		}
		std::cout << "Using option: " << optimalId.name << std::endl;
	}

	bool forceTune( int iter )
	{
		double bestPerformance = 0;
		bool tuneSuccessful = false;

		std::vector<RuntimeChooserOption>* ptrToOptions = getOptions();

		for( vector<RuntimeChooserOption>::iterator it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
		{
			try
			{
				std::cout << "warmup option " << it->name << ": ";
				measure( it->id, iter );
			}
			catch(InvalidKernelSetupException& e)
			{
				std::cout << "invalid configuration" << std::endl;
				continue;
			}
			std::cout << "done." << std::endl;
			break;
		}

		for( vector<RuntimeChooserOption>::iterator it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
		{
			preTuneAction();

			std::cout << "Tuning option " << it->name << ": " << std::flush;
			double performance;
			try
			{
				performance = 1./measure( it->id, iter );
				tuneSuccessful = true;
			}
			catch(InvalidKernelSetupException& e)
			{
				std::cout << "invalid configuration" << std::endl;
				continue;
			}
			CUDA_LAST_ERROR("TunableObject::measure()" );

			std::cout << performance << std::endl;

			if( performance > bestPerformance )
			{
				bestPerformance = performance;
				optimalId = *it;
			}
		}

		return tuneSuccessful;
	}

protected:
	RuntimeChooserOption optimalId;

private:
	Chronotimer timer;
};

}

#endif /* MEASURABLEOBJECT_H_ */
