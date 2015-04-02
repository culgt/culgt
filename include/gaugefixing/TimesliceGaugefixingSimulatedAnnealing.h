/**
 *  Created on: Apr 29, 2014
 *      Author: vogt
 *
 *
 * TODO this should be used by Coulomb gauge and by spatial maximal center gauge
 *
 * GaugeType = DirectMaximalCenterGaugeType<MCG_SPATIAL>
 * GaugeType = LandauCoulombGaugeType<COULOMB>
 */

#ifndef TIMESLICEGAUGEFIXINGSIMULATEDANNEALING_H_
#define TIMESLICEGAUGEFIXINGSIMULATEDANNEALING_H_

#include "CoulombGaugeTunableObject.h"
#include "algorithms/OrUpdate.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"
#include "TimesliceGaugeFixingOverrelaxationKernel.h"
#include "TimesliceGaugeFixingSimulatedAnnealingKernel.h"
#include "GaugefixingLaunchBounds.h"

#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"

namespace culgt
{

template<typename PatternType, typename LocalLinkType, typename GaugeType> class TimesliceGaugeFixingSimulatedAnnealing: public CoulombGaugeTunableObject<GlobalLink<PatternType,true>,LocalLinkType>
{
public:
	typedef GlobalLink<PatternType, true > GlobalLinkType;
	typedef CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef  GlobalLink<PatternType, UseTexture::value > GlobalLinkTypeInStep;
			COPY_GLOBALLINKTYPE( GlobalLinkTypeInStep, GlobalLinkTypeInStep2, 1 );

			KernelSetup<GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE::Ndim> setupSplit( object->super::dimTimeslice, true, GFLaunchBounds::SitesPerBlock );

			cudaFuncSetCacheConfig( TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

#if CUDART_VERSION >= 6050
			int numBlocks;
			size_t sharedMemorySize = 0; // TODO
			cudaError_t err = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &numBlocks, TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, GFLaunchBounds::SitesPerBlock*ThreadsPerSite::value, sharedMemorySize );
			if( numBlocks == 0 || err != cudaSuccess ) throw InvalidKernelSetupException();

			TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );

			CUDA_LAST_ERROR( "TimesliceGaugeFixingSimulatedAnnealing" );
#else
			TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );

			cudaDeviceSynchronize();
			cudaError_t error = cudaGetLastError();
			if( error == cudaErrorInvalidConfiguration || error == cudaErrorInvalidValue || error == cudaErrorLaunchOutOfResources ) // ugly...
			{
				// we encountered a wrong kernel setup (this is possible in autotune) exclude from auotune...
				throw InvalidKernelSetupException();
			}
			else
			{
				CUDA_ERROR( error, "TimesliceGaugeFixingSimulatedAnnealing" );
			}
#endif
		}
	};


	typedef TimesliceGaugeFixingSimulatedAnnealing<PatternType, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4 > threadsPerSiteSequence;
//	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	TimesliceGaugeFixingSimulatedAnnealing( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice, long seed, float temperature  ) : CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType>( Ut, UtDown, dimTimeslice, seed ), temperature(temperature)
	{
		Chooser::object = this;
	}

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor, bool UseTexture> inline void step()
	{
		Step<GaugefixingLaunchBounds<SitesPerBlock,MinBlocksPerMultiprocessor>, boost::mpl::int_<ThreadsPerSite>, boost::mpl::bool_<UseTexture> >::exec( this );
	}

	std::vector<RuntimeChooserOption>* getOptions()
	{
		return &Chooser::options;
	}

	void runNew( size_t id )
	{
		runner.run( id );
	}

	void run()
	{
		run( super::optimalId.id );
	}

	void run( size_t id )
	{
		runNew( id );
	}

	void setTemperature(float temperature)
	{
		this->temperature = temperature;
	}

private:
	float temperature;
};


}


#endif
