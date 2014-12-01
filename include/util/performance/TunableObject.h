/**
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef TUNABLEOBJECT_H_
#define TUNABLEOBJECT_H_

#include <boost/config.hpp>
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
#define CULGT_USE_CXX11_AUTOTUNE
#endif

#if __cplusplus < 201103L
#undef CULGT_USE_CXX11_AUTOTUNE
#endif

#ifdef CULGT_USE_CXX11_AUTOTUNE
#warning USING CX11 AUTOTUNE
#else
#warning USING OLD AUTOTUNE
#endif

#include "../../cuLGT1legacy/Chronotimer.h"
#include "AutotuneManager.h"

#include "cuda.h"
#include <vector>

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
//	virtual void runNew( size_t id ) = 0;
	virtual void preTuneAction() = 0;

#ifdef CULGT_USE_CXX11_AUTOTUNE
	virtual std::vector<size_t>* getOptions() = 0;
#endif

	virtual string getClassName()
	{
		return typeid( *this ).name();
	}

#ifdef CULGT_USE_CXX11_AUTOTUNE
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
#else
	double measure( int id, int iter = 100 )
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
#endif


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

#ifdef CULGT_USE_CXX11_AUTOTUNE
	void forceTune( int iter )
	{
		double bestPerformance = 0;

		std::vector<size_t>* ptrToOptions = getOptions();

		std::cout << ptrToOptions->size() << std::endl;

		for( auto it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
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

		for( auto it = ptrToOptions->begin(); it != ptrToOptions->end(); ++it )
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
#else
	void forceTune( int iter )
	{
		double bestPerformance = 0;
		int id = 0;

		// warmup

		int warmupId = 0;
		while( true )
		{
			try
			{
				measure( warmupId, iter );
			}
			catch(InvalidKernelSetupException& e)
			{
				warmupId++;
				continue;
			}
			catch(LastElementReachedException& e )
			{
				break;
			}
			std::cout << "warmup done..." << std::endl;
			break;
		}

		while( true )
		{
			try
			{
				preTuneAction();

				double performance;
				try
				{
					performance = 1./measure( id, iter );
				}
				catch(InvalidKernelSetupException& e)
				{
					std::cout << "Tuning option " << id << ": invalid configuration" << std::endl;
					id++;
					continue;
				}
				CUDA_LAST_ERROR("TunableObject::measure()" );

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
	}
#endif

protected:
	size_t optimalId;

private:
	Chronotimer timer;
};

}

#endif /* MEASURABLEOBJECT_H_ */
