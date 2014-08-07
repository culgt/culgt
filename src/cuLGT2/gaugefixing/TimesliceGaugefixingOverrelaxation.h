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

namespace culgt
{

namespace TimesliceGaugeFixingOverrelaxationKernel
{
	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter)
	{
		typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
		Algorithm orupdate( orParameter );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}


template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType> class TimesliceGaugeFixingOverrelaxation: public CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType>
{
public:
	typedef CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	TimesliceGaugeFixingOverrelaxation( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice, float orParameter, long seed ) : CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType>( Ut, UtDown, dimTimeslice, seed ), orParameter(orParameter)
	{
	}

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> inline void step()
	{
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( super::dimTimeslice, true, SitesPerBlock );

		cudaFuncSetCacheConfig( TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

		TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::Ut, *super::UtDown, super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dimTimeslice ), false, orParameter );
		TimesliceGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::Ut, *super::UtDown, super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dimTimeslice ), true, orParameter );
	}

	void run( int id = -1 )
	{
		if( id == -1 ) id = super::optimalId;

		switch( id )
		{
			case 0:
				step<8,32,1>();
				break;
			case 1:
				step<8,32,2>();
				break;
			case 2:
				step<8,32,3>();
				break;
			case 3:
				step<8,32,4>();
				break;
			case 4:
				step<8,32,5>();
				break;
			case 5:
				step<8,32,6>();
				break;
			case 6:
				step<8,64,1>();
				break;
			case 7:
				step<8,64,2>();
				break;
			case 8:
				step<8,64,3>();
				break;
			case 9:
				step<8,128,1>();
				break;
			case 10:
				step<4,32,1>();
				break;
			case 11:
				step<4,32,2>();
				break;
			case 12:
				step<4,32,3>();
				break;
			case 13:
				step<4,32,4>();
				break;
			case 14:
				step<4,32,5>();
				break;
			case 15:
				step<4,32,6>();
				break;
			case 16:
				step<4,64,1>();
				break;
			case 17:
				step<4,64,2>();
				break;
			case 18:
				step<4,64,3>();
				break;
			case 19:
				step<4,128,1>();
				break;
			default:
				throw LastElementReachedException();
		}
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
