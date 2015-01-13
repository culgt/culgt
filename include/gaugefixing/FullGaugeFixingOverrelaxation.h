/**
 * We use partial template specialization for the use GPUPatternParityPriority or GPUPatternTimesliceParityPriority.
 * The first choice is the natural one for Landau-like gauges, i.e. gauges that include all links. But the second choice
 * might be useful if you also deal with Coulomb-like gauges that operate only in timeslices. Then there is no need to transform
 * the configurations. Performance should only degrade in the latter case if we cannot fully utilize the GPU with timeslice operations.
 * TODO in the current implemtation of the timeslice case we bind/unbind textures in every iteration. That might be a performance obstacle.
 *
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef FULLGAUGEFIXINGOVERRELAXATION_H_
#define FULLGAUGEFIXINGOVERRELAXATION_H_

#include "LandauGaugeTunableObject.h"
#include "algorithms/OrUpdate.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"
#include "FullGaugeFixingOverrelaxationKernel.h"
#include "TimesliceGaugeFixingOverrelaxationKernel.h"
#include "../lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "GaugefixingLaunchBounds.h"

#include <iostream>

#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>
#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"

using boost::mpl::placeholders::_;

namespace culgt
{

template<typename PatternType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingOverrelaxation: public LandauGaugeTunableObject<GlobalLink<PatternType,true>,LocalLinkType>
{
public:
	typedef GlobalLink<PatternType,true> GlobalLinkType;
	typedef LandauGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	void run( int id = -1 )
	{
		assert(false);
	};
};

// disassemble the GlobalLink to specialize for the PatternType: here specialized for the GPUPatternParityPriority
template<typename SiteType, typename ParamType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingOverrelaxation<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType>: public LandauGaugeTunableObject<GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, true>,LocalLinkType>
{
public:
	typedef GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, true > GlobalLinkType;
	typedef LandauGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename ParamType::TYPE T;
	typedef typename ParamType::REALTYPE REALT;


	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, UseTexture::value > GlobalLinkTypeInStep;
			KernelSetup<SiteType::Ndim> setupSplit( object->super::dim, true, GFLaunchBounds::SitesPerBlock );

			cudaFuncSetCacheConfig( FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

			GlobalLinkTypeInStep::bindTexture( *object->super::U, GlobalLinkTypeInStep::getArraySize( object->super::dim ) );

			FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), false, object->orParameter );
			FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), true, object->orParameter );

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
		}
	};

	typedef FullGaugeFixingOverrelaxation<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<SiteType::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
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

	void run( size_t id  )
	{
		runNew( id );
	}

	void setOrParameter(float orParameter)
	{
		this->orParameter = orParameter;
	}

private:
	float orParameter;
};

template<typename SiteType, typename ParamType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingOverrelaxation<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType>: public LandauGaugeTunableObject<GlobalLink<GPUPatternTimesliceParityPriority<SiteType,ParamType>, true>,LocalLinkType>
{
public:
	typedef GlobalLink<GPUPatternTimesliceParityPriority<SiteType,ParamType>, true> GlobalLinkType;

	typedef SiteIndex<SiteType::Ndim, FULL_SPLIT> SiteTypeTimeslice;

	typedef LandauGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename ParamType::TYPE T;
	typedef typename ParamType::REALTYPE REALT;

	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef GlobalLink<GPUPatternParityPriority<SiteTypeTimeslice,ParamType>, true> GlobalLinkTypeTimeslice;
			COPY_GLOBALLINKTYPE( GlobalLinkTypeTimeslice, GlobalLinkTypeTimeslice2, 1 );

			KernelSetup<SiteType::Ndim> setupSplit( object->super::dim.getDimensionTimeslice(), true, GFLaunchBounds::SitesPerBlock );

			cudaFuncSetCacheConfig( TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

			int timesliceArraySize = GlobalLinkType::getArraySize( object->super::dim.getDimensionTimeslice() );
//			int t = 0;
			for( int t = 0; t < object->super::dim.getDimension(0); t++ )
			{
				int tdown = (t==0)?(object->super::dim.getDimension(0)-1):(t-1);

				// TODO this binding is most probably not a good idea... Better to use one bound and write a new kernel that deals with a full array by passing the timeslice...

				GlobalLinkTypeTimeslice::unbindTexture();
				GlobalLinkTypeTimeslice::bindTexture( &((*object->super::U)[t*timesliceArraySize]), timesliceArraySize );
				GlobalLinkTypeTimeslice2::unbindTexture();
				GlobalLinkTypeTimeslice2::bindTexture( &((*object->super::U)[tdown*timesliceArraySize]), timesliceArraySize);

				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), false, object->orParameter );
				TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), true, object->orParameter );
			}


			cudaDeviceSynchronize();
			cudaError_t error = cudaGetLastError();
			if( error != cudaSuccess )
			if( error == cudaErrorInvalidConfiguration || error == cudaErrorInvalidValue || error == cudaErrorLaunchOutOfResources )
			{
				// we encountered a wrong kernel setup (this is possible in autotune) exclude from auotune...
				throw InvalidKernelSetupException();
			}
			else
			{
				CUDA_ERROR( error, "FullGaugeFixingOverrelaxation" );
			}
		}
	};


	typedef FullGaugeFixingOverrelaxation<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<SiteType::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
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
		run( super::optimalId );
	}

	void run( size_t id  )
	{
		runNew( id );
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
