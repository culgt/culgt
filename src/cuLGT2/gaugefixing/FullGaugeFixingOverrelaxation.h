/**
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

namespace culgt
{

namespace FullGaugeFixingOverrelaxationKernel
{
	template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelOrStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE orParameter )
	{
		typedef OrUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
		Algorithm orupdate( orParameter );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, GaugeType, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( orupdate );
			gaugefixing.applyAlgorithm( U, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}


template<typename GlobalLinkType, typename LocalLinkType, typename GaugeType> class FullGaugeFixingOverrelaxation: public LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>
{
public:
	typedef LandauGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	FullGaugeFixingOverrelaxation( T** U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, float orParameter, long seed ) : LandauGaugeTunableObject<GlobalLinkType,LocalLinkType>( U, dim, seed ), orParameter(orParameter)
	{
	}

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> inline void step()
	{
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( super::dim, true, SitesPerBlock );

		cudaFuncSetCacheConfig( FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

		FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::U, super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dim ), false, orParameter );
		FullGaugeFixingOverrelaxationKernel::kernelOrStep<GlobalLinkType,LocalLinkType,GaugeType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,GaugeType::SharedArraySize*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::U, super::dim.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dim ), true, orParameter );
		CUDA_LAST_ERROR("FullGaugeFixingOverrelaxation" );
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



			case 20:
				step<4,32,7>();
				break;
			case 21:
				step<4,32,8>();
				break;

			case 22:
				step<4,32,9>();
				break;
			case 23:
				step<4,32,10>();
				break;
			case 24:
				step<4,32,11>();
				break;
			case 25:
				step<4,32,12>();
				break;
			case 26:
				step<4,64,4>();
				break;
			case 27:
				step<4,64,5>();
				break;
			case 28:
				step<4,64,6>();
				break;
			case 29:
				step<4,128,2>();
				break;
			case 30:
				step<4,128,3>();
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
