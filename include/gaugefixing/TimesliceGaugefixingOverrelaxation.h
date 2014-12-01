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

#include "CoulombGaugeTunableObject.h"
#include "algorithms/OrUpdate.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"
#include "TimesliceGaugeFixingOverrelaxationKernel.h"
#include "GaugefixingLaunchBounds.h"

#ifdef CULGT_USE_CXX11_AUTOTUNE
#include "../util/template_instantiation/SequenceRunner.h"
#include "../util/template_instantiation/RuntimeChooser.h"
#endif

namespace culgt
{

template<typename PatternType, typename LocalLinkType, typename GaugeType> class TimesliceGaugeFixingOverrelaxation: public CoulombGaugeTunableObject<GlobalLink<PatternType,true>,LocalLinkType>
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

			cudaFuncSetCacheConfig( TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

			TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), false, object->orParameter );
			TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkTypeInStep,GlobalLinkTypeInStep2,LocalLinkType,GaugeType,ThreadsPerSite::value,GFLaunchBounds::SitesPerBlock,GFLaunchBounds::MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite::value,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *object->super::Ut, *object->super::UtDown, object->super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( object->super::dimTimeslice ), true, object->orParameter );

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

	typedef TimesliceGaugeFixingOverrelaxation<PatternType, LocalLinkType, GaugeType> thisClass;
	typedef RuntimeChooser<thisClass, Step<_,_,_> > Chooser;
//	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<128,1>> launchBoundsSequence;
	typedef mpl::vector< GaugefixingLaunchBounds<32,1>, GaugefixingLaunchBounds<32,2>, GaugefixingLaunchBounds<32,3>, GaugefixingLaunchBounds<32,4>, GaugefixingLaunchBounds<32,5>, GaugefixingLaunchBounds<32,6>, GaugefixingLaunchBounds<32,7>, GaugefixingLaunchBounds<32,8>, GaugefixingLaunchBounds<64,1>, GaugefixingLaunchBounds<64,2>, GaugefixingLaunchBounds<64,3>, GaugefixingLaunchBounds<64,4>, GaugefixingLaunchBounds<64,5>, GaugefixingLaunchBounds<64,6>, GaugefixingLaunchBounds<128,1>, GaugefixingLaunchBounds<128,2>, GaugefixingLaunchBounds<128,3> > launchBoundsSequence;
	typedef mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
	typedef mpl::vector_c< int, 0, 1 > useTextureSequence;
	SequenceRunnerFrontend<Chooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> runner;
#endif

	TimesliceGaugeFixingOverrelaxation( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice, float orParameter, long seed ) : CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType>( Ut, UtDown, dimTimeslice, seed ), orParameter(orParameter)
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

	void run( size_t id )
	{
#ifdef CULGT_USE_CXX11_AUTOTUNE
		runNew( id );
#else
		switch( id )
		{
			case 0:
				step<8,32,1,false>();
				break;
			case 1:
				step<8,32,2,false>();
				break;
			case 2:
				step<8,32,3,false>();
				break;
			case 3:
				step<8,32,4,false>();
				break;
			case 4:
				step<8,32,5,false>();
				break;
			case 5:
				step<8,32,6,false>();
				break;
			case 6:
				step<8,64,1,false>();
				break;
			case 7:
				step<8,64,2,false>();
				break;
			case 8:
				step<8,64,3,false>();
				break;
			case 9:
				step<8,128,1,false>();
				break;
			case 10:
				step<4,32,1,false>();
				break;
			case 11:
				step<4,32,2,false>();
				break;
			case 12:
				step<4,32,3,false>();
				break;
			case 13:
				step<4,32,4,false>();
				break;
			case 14:
				step<4,32,5,false>();
				break;
			case 15:
				step<4,32,6,false>();
				break;
			case 16:
				step<4,64,1,false>();
				break;
			case 17:
				step<4,64,2,false>();
				break;
			case 18:
				step<4,64,3,false>();
				break;
			case 19:
				step<4,128,1,false>();
				break;

			case 20:
				step<8,32,1,true>();
				break;
			case 21:
				step<8,32,2,true>();
				break;
			case 22:
				step<8,32,3,true>();
				break;
			case 23:
				step<8,32,4,true>();
				break;
			case 24:
				step<8,32,5,true>();
				break;
			case 25:
				step<8,32,6,true>();
				break;
			case 26:
				step<8,64,1,true>();
				break;
			case 27:
				step<8,64,2,true>();
				break;
			case 28:
				step<8,64,3,true>();
				break;
			case 29:
				step<8,128,1,true>();
				break;
			case 30:
				step<4,32,1,true>();
				break;
			case 31:
				step<4,32,2,true>();
				break;
			case 32:
				step<4,32,3,true>();
				break;
			case 33:
				step<4,32,4,true>();
				break;
			case 34:
				step<4,32,5,true>();
				break;
			case 35:
				step<4,32,6,true>();
				break;
			case 36:
				step<4,64,1,true>();
				break;
			case 37:
				step<4,64,2,true>();
				break;
			case 38:
				step<4,64,3,true>();
				break;
			case 39:
				step<4,128,1,true>();
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
