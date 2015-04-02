/**
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXING4THREADS_H_
#define GAUGEFIXING4THREADS_H_
#include "../lattice/LocalLink.h"
#include "../lattice/parameterization_types/SU2Vector4.h"
#include "../lattice/SubgroupIterator.h"
#include "../lattice/GlobalLink.h"

using culgt::LocalLink;
using culgt::GlobalLink;
using culgt::SU2Vector4;

#include <boost/mpl/assert.hpp>

template<typename Algorithm, typename GaugeType, typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock=32> class GaugeFixing4Threads
{
public:
	__device__ inline GaugeFixing4Threads( Algorithm algorithm ) : algorithm(algorithm), id(threadIdx.x % SitesPerBlock), mu(threadIdx.x / SitesPerBlock)
	{
	}

	__device__ inline void applyAlgorithm( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t* nn, lat_index_t latticeSize, bool parity )
	{
		// we do not need the lattice extents in all directions in this kernel: we can save registers by only giving the total size
		typename GlobalLinkType::PATTERNTYPE::SITETYPE site(latticeSize, nn );

		lat_index_t index = blockIdx.x * blockDim.x/threadsPerSite + id;

		if( index >= latticeSize ) return;

		if( parity == 1 ) index += site.getSize()/2;

		site.setIndex( index );

		GlobalLinkType globalLinkUp( U, site, mu );
		localLinkUp = globalLinkUp;

		site.setNeighbour( mu, false );
		GlobalLinkType globalLinkDown( U, site, mu );
		localLinkDown = globalLinkDown;

		SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

		globalLinkUp = localLinkUp;
		globalLinkDown = localLinkDown;
	}

	template<typename GlobalLinkType2> __device__ inline void applyAlgorithmTimeslice( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Uup, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Udown , lat_index_t* nn, lat_index_t latticeSize, bool parity )
	{
		// we do not need the lattice extents in all directions in this kernel: we can save registers by only giving the total size

		typename GlobalLinkType::PATTERNTYPE::SITETYPE site(latticeSize, nn );
		lat_index_t index = blockIdx.x * blockDim.x/threadsPerSite + id;


		if( parity == 1 ) index += site.getSize()/2;

		site.setIndex( index );

		GlobalLinkType globalLinkUp( Uup, site, mu );
		localLinkUp = globalLinkUp;

		if( mu != 0 )
		{
			site.setNeighbour( mu, false );

			GlobalLinkType globalLinkDown( Uup, site, mu );
			localLinkDown = globalLinkDown;

			SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

			globalLinkUp = localLinkUp;
			globalLinkDown = localLinkDown;
		}
		else
		{
//			COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );
			GlobalLinkType2 globalLinkDown( Udown, site, mu );
			localLinkDown = globalLinkDown;

			SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

			globalLinkUp = localLinkUp;
			globalLinkDown = localLinkDown;
		}

	}

	__device__ inline void subgroupStep( lat_group_index_t iSub, lat_group_index_t jSub )
	{
		extern __shared__ typename LocalLinkType::PARAMTYPE::REALTYPE shA[]; // define size in kernel call (size needs to be 4*NSB)!
		initializeSharedMemory( shA, GaugeType::SharedArraySize );

		__syncthreads();
		LocalLink<SU2Vector4<typename LocalLinkType::PARAMTYPE::REALTYPE> > quaternionUp = localLinkUp.getSU2Subgroup( iSub, jSub );
		LocalLink<SU2Vector4<typename LocalLinkType::PARAMTYPE::REALTYPE> > quaternionDown = localLinkDown.getSU2Subgroup( iSub, jSub );

		typename Real4<typename LocalLinkType::PARAMTYPE::REALTYPE>::VECTORTYPE& qUp = quaternionUp[0];
		typename Real4<typename LocalLinkType::PARAMTYPE::REALTYPE>::VECTORTYPE& qDown = quaternionDown[0];

		GaugeType::gatherInfo( shA, qUp, qDown, id, mu, SitesPerBlock );

		// calc update
		__syncthreads();
		if( mu == 0 )
		{
			GaugeType::prepareUpdate( shA, id, SitesPerBlock );
			algorithm.calculateUpdate( shA, id, SitesPerBlock );
		}
		__syncthreads();

		GaugeType::collectUpdate( shA, qUp, id, mu, 0, SitesPerBlock );
		__syncthreads(); // this is necessary otherwise the next subgroup already overwrites shA!

		quaternionUp.set( 0, qUp );


		localLinkUp.leftSubgroupMult( quaternionUp, iSub, jSub );
		quaternionUp.hermitian();
		localLinkDown.rightSubgroupMult( quaternionUp, iSub, jSub );
	}
private:
	static const int threadsPerSite = 4;
	const short id;
	const short mu;
	Algorithm algorithm;
	LocalLinkType localLinkUp;
	LocalLinkType localLinkDown;

	__device__ inline void initializeSharedMemory( typename LocalLinkType::PARAMTYPE::REALTYPE* shA, int sharedArraySize )
	{
		if( mu == 0 )
		{
			for( int i = 0; i < sharedArraySize; i++ )
			{
				shA[id+i*SitesPerBlock] = 0;
			}
		}
		__syncthreads();
	}
};

#endif /* GAUGEFIXINGSUBGROUPSTEP_H_ */
