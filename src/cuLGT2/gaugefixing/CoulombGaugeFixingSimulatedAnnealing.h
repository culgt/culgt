/**
 *
 *  Created on: Apr 29, 2014
 *      Author: vogt
 */

#ifndef COULOMBGAUGEFIXINGSIMULATEDANNEALING_H_
#define COULOMBGAUGEFIXINGSIMULATEDANNEALING_H_

#include "CoulombGaugeTunableObject.h"
#include "algorithms/SaUpdate.h"
#include "GaugeFixing8Threads.h"
#include "GaugeFixing4Threads.h"
#include "../util/rng/PhiloxWrapper.h"

namespace culgt
{

namespace CoulombGaugeFixingSimulatedAnnealingKernel
{
	template<typename GlobalLinkType, typename LocalLinkType, int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> __global__ __launch_bounds__(ThreadsPerSite*SitesPerBlock,MinBlocksPerMultiprocessor) void kernelSaStep( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Ut, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* UtDown, lat_index_t latticeSize, lat_index_t* nn, bool parity, typename LocalLinkType::PARAMTYPE::REALTYPE temperature, int rngSeed, int rngCounter)
	{
		PhiloxWrapper<typename LocalLinkType::PARAMTYPE::REALTYPE> rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
		typedef SaUpdate<typename LocalLinkType::PARAMTYPE::REALTYPE> Algorithm;
		Algorithm saupdate( temperature, rng );

		if( ThreadsPerSite == 8 )
		{
			GaugeFixing8Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
			gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
		}
		else if( ThreadsPerSite == 4 )
		{
			GaugeFixing4Threads<Algorithm, LandauCoulombGaugeType<COULOMB>, GlobalLinkType, LocalLinkType, SitesPerBlock> gaugefixing( saupdate );
			gaugefixing.applyAlgorithmTimeslice( Ut, UtDown, nn, latticeSize, parity );
		}
		else
		{
			assert( false );
		}
	}
}


template<typename GlobalLinkType, typename LocalLinkType> class CoulombGaugeFixingSimulatedAnnealing: public CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType>
{
public:
	typedef CoulombGaugeTunableObject<GlobalLinkType,LocalLinkType> super;

	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	inline CoulombGaugeFixingSimulatedAnnealing( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dimTimeslice, long seed, float temperature ) : super(Ut,UtDown,dimTimeslice,seed), temperature(temperature)
	{
	}

	template<int ThreadsPerSite, int SitesPerBlock, int MinBlocksPerMultiprocessor> inline void step()
	{
		KernelSetup<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> setupSplit( super::dimTimeslice, true, SitesPerBlock );

		cudaFuncSetCacheConfig( CoulombGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkType,LocalLinkType,ThreadsPerSite, SitesPerBlock,MinBlocksPerMultiprocessor>, cudaFuncCachePreferL1 );

		CoulombGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkType,LocalLinkType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::Ut, *super::UtDown, super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dimTimeslice ), false, temperature, super::seed, PhiloxWrapper<REALT>::getNextCounter() );
		CoulombGaugeFixingSimulatedAnnealingKernel::kernelSaStep<GlobalLinkType,LocalLinkType,ThreadsPerSite,SitesPerBlock,MinBlocksPerMultiprocessor><<<setupSplit.getGridSize(),setupSplit.getBlockSize()*ThreadsPerSite,4*setupSplit.getBlockSize()*sizeof(REALT)>>>( *super::Ut, *super::UtDown, super::dimTimeslice.getSize(), SiteNeighbourTableManager<typename GlobalLinkType::PATTERNTYPE::SITETYPE>::getDevicePointer( super::dimTimeslice ), true, temperature, super::seed, PhiloxWrapper<REALT>::getNextCounter() );
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

	void setTemperature(float temperature)
	{
		this->temperature = temperature;
	}

private:
	float temperature;
};


}


#endif /* COULOMBGAUGEFIXINGOVERRELAXATION_H_ */
