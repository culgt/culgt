/**
 * We use partial template specialization for choosing between GPUPatternParityPriority or GPUPatternTimesliceParityPriority.
 * The first choice is the natural one for Landau-like gauges, i.e. gauges that include all links. But the second choice
 * might be useful if you also deal with Coulomb-like gauges that operate only in timeslices. Then there is no need to transform
 * the configurations. Performance should only degrade in the latter case if we cannot fully utilize the GPU with timeslice operations.
 */

#ifndef FULLGAUGEFIXINGSIMULATEDANNEALING_H_
#define FULLGAUGEFIXINGSIMULATEDANNEALING_H_

#include "FullGaugeTunableObject.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"
#include "../lattice/configuration_patterns/GPUPatternTimesliceParityPriority.h"
#include "GaugefixingLaunchBounds.h"
#include "FullGaugeFixingSimulatedAnnealingKernel.h"
#include "TimesliceGaugeFixingSimulatedAnnealingKernel.h"

#include <iostream>

#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>
#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"

using boost::mpl::placeholders::_;

namespace culgt
{

/**
 * General class, needs specialization for different Patterns
 */
template<typename PatternType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingSimulatedAnnealing: public FullGaugeTunableObject<GlobalLink<PatternType,true>,LocalLinkType>
{
public:
	typedef GlobalLink<PatternType,true> GlobalLinkType;
	typedef FullGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	void run( int id = -1 )
	{
		assert(false);
	};
};

/**
 * Specialization for PatternType == GPUPatternParityPriority (i.e. the natural choice for Landau-like gauges)
 */
template<typename SiteType, typename ParamType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingSimulatedAnnealing<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType>: public FullGaugeTunableObject<GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, true>,LocalLinkType>
{
public:
	typedef GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, true > GlobalLinkType;
	typedef FullGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename ParamType::TYPE T;
	typedef typename ParamType::REALTYPE REALT;

	/**
	 * Perform a gaugefixing step with the given kernel settings.
	 */
	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef GlobalLink<GPUPatternParityPriority<SiteType,ParamType>, UseTexture::value > GlobalLinkTypeInStep;
			KernelSetup<SiteType::NDIM> setupSplit( object->super::dim, true, GFLaunchBounds::SitesPerBlock );

			cudaFuncSetCacheConfig( FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

			GlobalLinkTypeInStep::bindTexture( *object->super::U, GlobalLinkTypeInStep::getArraySize( object->super::dim ) );

#if CUDART_VERSION >= 6050
			int numBlocks;
			size_t sharedMemorySize = GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT);
			cudaError_t err = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &numBlocks, FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, GFLaunchBounds::SitesPerBlock*ThreadsPerSite::value, sharedMemorySize );
			if( numBlocks == 0 || err != cudaSuccess ) throw InvalidKernelSetupException();

			FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			CUDA_LAST_ERROR( "FullGaugeFixingSimulatedAnnealing" );
#else
			FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			FullGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeInStep,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::U, object->super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkTypeInStep::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dim ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );

			cudaDeviceSynchronize();
			cudaError_t error = cudaGetLastError();
			if( error == cudaErrorInvalidConfiguration || error == cudaErrorInvalidValue || error == cudaErrorLaunchOutOfResources ) // ugly...
			{
				// we encountered a wrong kernel setup (this is possible in autotune) exclude from auotune...
				throw InvalidKernelSetupException();
			}
			else
			{
				CUDA_ERROR( error, "FullGaugeFixingSimulatedAnnealing" );
			}
#endif
		}
	};

	typedef FullGaugeFixingSimulatedAnnealing<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;

#ifdef CULGT_NO_AUTOTUNE
#warning "NOT using autotune for FullGaugeFixingSimulatedAnnealing!"
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,6> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 1 > useTextureSequence;

	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingSimulatedAnnealing( T** U, LatticeDimension<SiteType::NDIM> dim, long seed, float temperature ) : FullGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), temperature(temperature), runner( RUN_FIRST_CHOICE )
	{
		Chooser::object = this;
	}
#else
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;

	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingSimulatedAnnealing( T** U, LatticeDimension<SiteType::NDIM> dim, long seed, float temperature ) : FullGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), temperature(temperature)
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

	void run( size_t id  )
	{
		runner.run( id );
	}

	void setTemperature(float temperature)
	{
		this->temperature = temperature;
	}

private:
	float temperature;
};

/**
 * Specialization for PatternType == GPUPatternTimesliceParityPriority (i.e. the natural choice for Coulomb-like gauges)
 */
template<typename SiteType, typename ParamType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingSimulatedAnnealing<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType>: public FullGaugeTunableObject<GlobalLink<GPUPatternTimesliceParityPriority<SiteType,ParamType>, true>,LocalLinkType>
{
public:
	typedef GlobalLink<GPUPatternTimesliceParityPriority<SiteType,ParamType>, true> GlobalLinkType;
	typedef SiteIndex<SiteType::NDIM, FULL_SPLIT> SiteTypeTimeslice;
	typedef FullGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename ParamType::TYPE T;
	typedef typename ParamType::REALTYPE REALT;

	/**
	 * Perform a gaugefixing step with the given kernel settings.
	 */
	template<typename GFLaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Step
	{
		template<typename T> static void exec( T* object )
		{
			typedef GlobalLink<GPUPatternParityPriority<SiteTypeTimeslice,ParamType>, UseTexture::value> GlobalLinkTypeTimeslice;
			COPY_GLOBALLINKTYPE( GlobalLinkTypeTimeslice, GlobalLinkTypeTimeslice2, 1 );

			KernelSetup<SiteType::NDIM> setupSplit( object->super::dim.getDimensionTimeslice(), true, GFLaunchBounds::SitesPerBlock );

			cudaFuncSetCacheConfig( TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

#if CUDART_VERSION >= 6050
			int numBlocks;
			size_t sharedMemorySize = GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT);
			cudaError_t err = cudaOccupancyMaxActiveBlocksPerMultiprocessor( &numBlocks, TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, GFLaunchBounds::SitesPerBlock*ThreadsPerSite::value, sharedMemorySize );
			if( numBlocks == 0 || err != cudaSuccess ) throw InvalidKernelSetupException();

			int timesliceArraySize = GlobalLinkType::getArraySize( object->super::dim.getDimensionTimeslice() );
			for( int t = 0; t < object->super::dim.getDimension(0); t++ )
			{
				int tdown = (t==0)?(object->super::dim.getDimension(0)-1):(t-1);

				// This binding in every step is most probably not a good idea... Better to use one bind to the full array and write a new kernel that deals with a full array by passing the timeslice...
				// Probably we won't use this choice very much. Until then we don't care...
				GlobalLinkTypeTimeslice::unbindTexture();
				GlobalLinkTypeTimeslice::bindTexture( &((*object->super::U)[t*timesliceArraySize]), timesliceArraySize );
				GlobalLinkTypeTimeslice2::unbindTexture();
				GlobalLinkTypeTimeslice2::bindTexture( &((*object->super::U)[tdown*timesliceArraySize]), timesliceArraySize);

				TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
				TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
			}
			CUDA_LAST_ERROR( "FullGaugeFixingSimulatedAnnealing-TimeslicePattern" );
#else
			int timesliceArraySize = GlobalLinkType::getArraySize( object->super::dim.getDimensionTimeslice() );
			for( int t = 0; t < object->super::dim.getDimension(0); t++ )
			{
				int tdown = (t==0)?(object->super::dim.getDimension(0)-1):(t-1);

				// This binding in every step is most probably not a good idea... Better to use one bind to the full array and write a new kernel that deals with a full array by passing the timeslice...
				// Probably we won't use this choice very much. Until then we don't care...
				GlobalLinkTypeTimeslice::unbindTexture();
				GlobalLinkTypeTimeslice::bindTexture( &((*object->super::U)[t*timesliceArraySize]), timesliceArraySize );
				GlobalLinkTypeTimeslice2::unbindTexture();
				GlobalLinkTypeTimeslice2::bindTexture( &((*object->super::U)[tdown*timesliceArraySize]), timesliceArraySize);

				TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), false, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
				TimesliceGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkTypeTimeslice,GlobalLinkTypeTimeslice2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( &((*object->super::U)[t*timesliceArraySize]), &((*object->super::U)[tdown*timesliceArraySize]), object->super::dim.getSizeTimeslice(), SiteNeighbourTableManager<SiteTypeTimeslice>::getDevicePointer( object->super::dim.getDimensionTimeslice() ), true, object->temperature, object->super::seed, PhiloxWrapper<REALT>::getNextCounter() );
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
				CUDA_ERROR( error, "FullGaugeFixingSimulatedAnnealing-TimeslicePattern" );
			}
#endif
		}
	};

	typedef FullGaugeFixingSimulatedAnnealing<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;

#ifdef CULGT_NO_AUTOTUNE
#warning "NOT using autotune for FullGaugeFixingSimulatedAnnealing-TimeslicePattern!"
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,6> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 1 > useTextureSequence;

	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingSimulatedAnnealing( T** U, LatticeDimension<SiteType::NDIM> dim, long seed, float temperature ) : FullGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), temperature(temperature), runner( RUN_FIRST_CHOICE )
	{
		Chooser::object = this;
	}
#else
	typedef boost::mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef boost::mpl::vector_c< int, 4 > threadsPerSiteSequence;
//	typedef boost::mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef boost::mpl::vector_c< int, 0, 1 > useTextureSequence;

	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingSimulatedAnnealing( T** U, LatticeDimension<SiteType::NDIM> dim, long seed, float temperature ) : FullGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), temperature(temperature)
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

	void run( size_t id  )
	{
		runner.run( id );
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
