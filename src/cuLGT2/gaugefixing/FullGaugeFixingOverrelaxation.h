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


#ifdef CULGT_USE_CXX11_AUTOTUNE

#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>
#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"

namespace mpl = boost::mpl;
using mpl::placeholders::_;
#endif

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

#ifdef CULGT_USE_CXX11_AUTOTUNE

	typedef FullGaugeFixingOverrelaxation<GPUPatternParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;
	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;
#endif

	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<SiteType::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
	{
#ifdef CULGT_USE_CXX11_AUTOTUNE
		Chooser::object = this;
#endif
	}

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor, bool UseTexture> inline void step()
	{
		Step<GaugefixingLaunchBounds<SitesPerBlock,MinBlocksPerMultiprocessor>, boost::mpl::int_<ThreadsPerSite>, boost::mpl::bool_<UseTexture> >::exec( this );
	}

#ifdef CULGT_USE_CXX11_AUTOTUNE
	std::vector<size_t>* getOptions()
	{
		return &Chooser::ids;
	}

	void runNew( size_t id )
	{
		runner.run( id );
	}
#endif

	void run()
	{
		run( super::optimalId );
	}

	void run( size_t id  )
	{
#ifdef CULGT_USE_CXX11_AUTOTUNE
		runNew( id );
#else
		switch( id )
		{
			case 0:
				step<8,32,1,true>();
				break;
			case 1:
				step<8,32,2,true>();
				break;
			case 2:
				step<8,32,3,true>();
				break;
			case 3:
				step<8,32,4,true>();
				break;
			case 4:
				step<8,32,5,true>();
				break;
			case 5:
				step<8,32,6,true>();
				break;
			case 6:
				step<8,64,1,true>();
				break;
			case 7:
				step<8,64,2,true>();
				break;
			case 8:
				step<8,64,3,true>();
				break;
			case 9:
				step<8,128,1,true>();
				break;
			case 10:
				step<4,32,1,true>();
				break;
			case 11:
				step<4,32,2,true>();
				break;
			case 12:
				step<4,32,3,true>();
				break;
			case 13:
				step<4,32,4,true>();
				break;
			case 14:
				step<4,32,5,true>();
				break;
			case 15:
				step<4,32,6,true>();
				break;
			case 16:
				step<4,64,1,true>();
				break;
			case 17:
				step<4,64,2,true>();
				break;
			case 18:
				step<4,64,3,true>();
				break;
			case 19:
				step<4,128,1,true>();
				break;



			case 20:
				step<4,32,7,true>();
				break;
			case 21:
				step<4,32,8,true>();
				break;

			case 22:
				step<4,32,9,true>();
				break;
			case 23:
				step<4,32,10,true>();
				break;
			case 24:
				step<4,32,11,true>();
				break;
			case 25:
				step<4,32,12,true>();
				break;
			case 26:
				step<4,64,4,true>();
				break;
			case 27:
				step<4,64,5,true>();
				break;
			case 28:
				step<4,64,6,true>();
				break;
			case 29:
				step<4,128,2,true>();
				break;
			case 30:
				step<4,128,3,true>();
				break;


			case 61:
				step<8,32,1,false>(); // <--- 61 here
				break;
			case 31:
				step<8,32,2,false>();
				break;
			case 32:
				step<8,32,3,false>();
				break;
			case 33:
				step<8,32,4,false>();
				break;
			case 34:
				step<8,32,5,false>();
				break;
			case 35:
				step<8,32,6,false>();
				break;
			case 36:
				step<8,64,1,false>();
				break;
			case 37:
				step<8,64,2,false>();
				break;
			case 38:
				step<8,64,3,false>();
				break;
			case 39:
				step<8,128,1,false>();
				break;
			case 40:
				step<4,32,1,false>();
				break;
			case 41:
				step<4,32,2,false>();
				break;
			case 42:
				step<4,32,3,false>();
				break;
			case 43:
				step<4,32,4,false>();
				break;
			case 44:
				step<4,32,5,false>();
				break;
			case 45:
				step<4,32,6,false>();
				break;
			case 46:
				step<4,64,1,false>();
				break;
			case 47:
				step<4,64,2,false>();
				break;
			case 48:
				step<4,64,3,false>();
				break;
			case 49:
				step<4,128,1,false>();
				break;



			case 50:
				step<4,32,7,false>();
				break;
			case 51:
				step<4,32,8,false>();
				break;

			case 52:
				step<4,32,9,false>();
				break;
			case 53:
				step<4,32,10,false>();
				break;
			case 54:
				step<4,32,11,false>();
				break;
			case 55:
				step<4,32,12,false>();
				break;
			case 56:
				step<4,64,4,false>();
				break;
			case 57:
				step<4,64,5,false>();
				break;
			case 58:
				step<4,64,6,false>();
				break;
			case 59:
				step<4,128,2,false>();
				break;
			case 60:
				step<4,128,3,false>();
				break;
			default:
				throw LastElementReachedException();
		}
#endif
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


#ifdef CULGT_USE_CXX11_AUTOTUNE

	typedef FullGaugeFixingOverrelaxation<GPUPatternTimesliceParityPriority<SiteType,ParamType>, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;
	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;

	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<SiteType::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
	{
		Chooser::object = this;
	}
#else
	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<SiteType::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
	{
	}
#endif

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor, bool UseTexture> inline void step()
	{
		Step<GaugefixingLaunchBounds<SitesPerBlock,MinBlocksPerMultiprocessor>, boost::mpl::int_<ThreadsPerSite>, boost::mpl::bool_<UseTexture> >::exec( this );
	}

#ifdef CULGT_USE_CXX11_AUTOTUNE
	std::vector<size_t>* getOptions()
	{
		return &Chooser::ids;
	}

	void runNew( size_t id )
	{
		runner.run( id );
	}
#endif

	void run()
	{
		run( super::optimalId );
	}

	void run( size_t id  )
	{
#ifdef CULGT_USE_CXX11_AUTOTUNE
		runNew( id );
#else
		switch( id )
		{
		case 0:
			step<8,32,1,true>();
			break;
		case 1:
			step<8,32,2,true>();
			break;
		case 2:
			step<8,32,3,true>();
			break;
		case 3:
			step<8,32,4,true>();
			break;
		case 4:
			step<8,32,5,true>();
			break;
		case 5:
			step<8,32,6,true>();
			break;
		case 6:
			step<8,64,1,true>();
			break;
		case 7:
			step<8,64,2,true>();
			break;
		case 8:
			step<8,64,3,true>();
			break;
		case 9:
			step<8,128,1,true>();
			break;
		case 10:
			step<4,32,1,true>();
			break;
		case 11:
			step<4,32,2,true>();
			break;
		case 12:
			step<4,32,3,true>();
			break;
		case 13:
			step<4,32,4,true>();
			break;
		case 14:
			step<4,32,5,true>();
			break;
		case 15:
			step<4,32,6,true>();
			break;
		case 16:
			step<4,64,1,true>();
			break;
		case 17:
			step<4,64,2,true>();
			break;
		case 18:
			step<4,64,3,true>();
			break;
		case 19:
			step<4,128,1,true>();
			break;



		case 20:
			step<4,32,7,true>();
			break;
		case 21:
			step<4,32,8,true>();
			break;

		case 22:
			step<4,32,9,true>();
			break;
		case 23:
			step<4,32,10,true>();
			break;
		case 24:
			step<4,32,11,true>();
			break;
		case 25:
			step<4,32,12,true>();
			break;
		case 26:
			step<4,64,4,true>();
			break;
		case 27:
			step<4,64,5,true>();
			break;
		case 28:
			step<4,64,6,true>();
			break;
		case 29:
			step<4,128,2,true>();
			break;
		case 30:
			step<4,128,3,true>();
			break;


		case 61:
			step<8,32,1,false>(); // <--- 61 here
			break;
		case 31:
			step<8,32,2,false>();
			break;
		case 32:
			step<8,32,3,false>();
			break;
		case 33:
			step<8,32,4,false>();
			break;
		case 34:
			step<8,32,5,false>();
			break;
		case 35:
			step<8,32,6,false>();
			break;
		case 36:
			step<8,64,1,false>();
			break;
		case 37:
			step<8,64,2,false>();
			break;
		case 38:
			step<8,64,3,false>();
			break;
		case 39:
			step<8,128,1,false>();
			break;
		case 40:
			step<4,32,1,false>();
			break;
		case 41:
			step<4,32,2,false>();
			break;
		case 42:
			step<4,32,3,false>();
			break;
		case 43:
			step<4,32,4,false>();
			break;
		case 44:
			step<4,32,5,false>();
			break;
		case 45:
			step<4,32,6,false>();
			break;
		case 46:
			step<4,64,1,false>();
			break;
		case 47:
			step<4,64,2,false>();
			break;
		case 48:
			step<4,64,3,false>();
			break;
		case 49:
			step<4,128,1,false>();
			break;



		case 50:
			step<4,32,7,false>();
			break;
		case 51:
			step<4,32,8,false>();
			break;

		case 52:
			step<4,32,9,false>();
			break;
		case 53:
			step<4,32,10,false>();
			break;
		case 54:
			step<4,32,11,false>();
			break;
		case 55:
			step<4,32,12,false>();
			break;
		case 56:
			step<4,64,4,false>();
			break;
		case 57:
			step<4,64,5,false>();
			break;
		case 58:
			step<4,64,6,false>();
			break;
		case 59:
			step<4,128,2,false>();
			break;
		case 60:
			step<4,128,3,false>();
			break;
		default:
			throw LastElementReachedException();
		}
#endif
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
