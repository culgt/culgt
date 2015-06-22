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

#ifndef TIMESLICEGAUGEFIXINGOVERRELAXATION_H_
#define TIMESLICEGAUGEFIXINGOVERRELAXATION_H_

#include "TimesliceGaugeTunableObject.h"
#include "algorithms/OrUpdate.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"
#include "TimesliceGaugeFixingOverrelaxationKernel.h"
#include "GaugefixingLaunchBounds.h"

#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"

namespace culgt
{

template<typename PatternType, typename LocalLinkType, typename GaugeType, bool DoMicro = false> class TimesliceGaugeFixingOverrelaxation: public TimesliceGaugeTunableObject<GlobalLink<PatternType,true>,LocalLinkType>
{
public:
	typedef GlobalLink<PatternType, true > GlobalLinkType;
	typedef TimesliceGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;

	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef  GlobalLink<PatternType, UseTexture::value > GlobalLinkTypeInStep;
			COPY_GLOBALLINKTYPE( GlobalLinkTypeInStep, GlobalLinkTypeInStep2, 1 );

			KernelSetup<GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE::NDIM> setupSplit( object->super::dimTimeslice, true, GFLaunchBounds::SitesPerBlock );

			if( DoMicro )
				cudaFuncSetCacheConfig( TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );
			else
				cudaFuncSetCacheConfig( TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

#if CUDART_VERSION >= 6050
			int numBlocks;
			size_t sharedMemorySize = 0; // TODO


			if( DoMicro )
			{
				cudaError_t err = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &numBlocks, TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, GFLaunchBounds::SitesPerBlock*ThreadsPerSite::value, sharedMemorySize );
				if( numBlocks == 0 || err != cudaSuccess ) throw InvalidKernelSetupException();
				TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false );
				CUDA_LAST_ERROR( "TimesliceGaugeFixingOverrelaxation(DoMicro) parity = 0" );
				TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true );
				CUDA_LAST_ERROR( "TimesliceGaugeFixingOverrelaxation(DoMicro) parity = 1" );
			}
			else
			{
				cudaError_t err = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &numBlocks, TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, GFLaunchBounds::SitesPerBlock*ThreadsPerSite::value, sharedMemorySize );
				if( numBlocks == 0 || err != cudaSuccess ) throw InvalidKernelSetupException();
				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false, object->orParameter );
				CUDA_LAST_ERROR( "TimesliceGaugeFixingOverrelaxation parity = 0" );
				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true, object->orParameter );
				CUDA_LAST_ERROR( "TimesliceGaugeFixingOverrelaxation parity = 1" );
			}


#else
			if( DoMicro )
			{
				TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false );
				TimesliceGaugeFixingOverrelaxationKernel::kernelMicroStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true );
			}
			else
			{
				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false, object->orParameter );
				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true, object->orParameter );
			}

			cudaDeviceSynchronize();
			cudaError_t error = cudaGetLastError();
			if( error == cudaErrorInvalidConfiguration || error == cudaErrorInvalidValue || error == cudaErrorLaunchOutOfResources ) // ugly...
			{
				// we encountered a wrong kernel setup (this is possible in autotune) exclude from auotune...
				throw InvalidKernelSetupException();
			}
			else
			{
				CUDA_ERROR( error, "FullGaugeFixingOverrelaxation" );
			}
#endif
		}
	};


	typedef TimesliceGaugeFixingOverrelaxation<PatternType, LocalLinkType, GaugeType, DoMicro> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;

	//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;

#ifdef CULGT_NO_AUTOTUNE
#warning "NOT using autotune for TimesliceGaugeFixingOverrelaxation!"
	typedef boost::mpl::vector< GaugefixingLaunchBounds<64,4> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;
	TimesliceGaugeFixingOverrelaxation( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dimTimeslice, long seed, float orParameter = 1.0  ) : TimesliceGaugeTunableObject<GlobalLinkType,LocalLinkType>( Ut, UtDown, dimTimeslice, seed ), orParameter(orParameter), runner( RUN_FIRST_CHOICE )
	{
		Chooser::object = this;
	}

#else
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	TimesliceGaugeFixingOverrelaxation( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dimTimeslice, long seed, float orParameter = 1.0  ) : TimesliceGaugeTunableObject<GlobalLinkType,LocalLinkType>( Ut, UtDown, dimTimeslice, seed ), orParameter(orParameter)
	{
		Chooser::object = this;
	}
#endif



	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor, bool UseTexture> inline void step()
	{
		Step<GaugefixingLaunchBounds<SitesPerBlock,MinBlocksPerMultiprocessor>, boost::mpl::int_<ThreadsPerSite>, boost::mpl::bool_<UseTexture> >::exec( this );
	}


	std::vector<RuntimeChooserOption>* getOptions()
	{
		return &Chooser::options;
	}

	void run()
	{
		run( super::optimalId.id );
	}

	void run( size_t id )
	{
		runner.run( id );
	}

	void setOrParameter(float orParameter)
	{
		this->orParameter = orParameter;
	}

private:
	float orParameter;
};


}


#endif
