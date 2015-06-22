/**
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXING8THREADS_H_
#define GAUGEFIXING8THREADS_H_
#include "../lattice/LocalLink.h"
#include "../lattice/parameterization_types/SU2Vector4.h"
#include "../lattice/GlobalLink.h"

using culgt::LocalLink;
using culgt::SU2Vector4;
using culgt::GlobalLink;


template<typename Algorithm, typename GaugeType, typename GlobalLinkType, typename LocalLinkType, int SitesPerBlock = 32> class GaugeFixing8Threads
{
public:
	__device__ inline GaugeFixing8Threads( Algorithm algorithm ) : algorithm(algorithm), updown( threadIdx.x / (4*SitesPerBlock) ), mu((threadIdx.x % (4*SitesPerBlock)) / SitesPerBlock), id((threadIdx.x % (4*SitesPerBlock)) % SitesPerBlock)
	{
	}


	__device__ inline void applyAlgorithm( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* U, lat_index_t* nn, lat_index_t latticeSize, bool parity )
	{
		// we do not need the lattice extents in all directions in this kernel: we can save registers by only giving the total size
		typename GlobalLinkType::PATTERNTYPE::SITETYPE site( latticeSize, nn );

		lat_index_t index = blockIdx.x * blockDim.x/threadsPerSite + id;

		if( parity == 1 ) index += site.getSize()/2;

		site.setIndex( index );

		if( updown == 1 )
		{
			site.setNeighbour( mu, false );
		}

		// load Link
		GlobalLinkType globalLink( U, site, mu );
		localLink = globalLink;

		SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

		// write back Link
		globalLink = localLink;
	}

	template<typename GlobalLinkType2> __device__ inline void applyAlgorithmTimeslice( typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Uup, typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE* Udown , lat_index_t* nn, lat_index_t latticeSize, bool parity )
	{
		typename GlobalLinkType::PATTERNTYPE::SITETYPE site( latticeSize, nn );

		lat_index_t index = blockIdx.x * blockDim.x/threadsPerSite + id;

		if( parity == 1 ) index += site.getSize()/2;

		site.setIndex( index );

		if( updown == 1 && mu > 0 )
		{
			site.setNeighbour( mu, false );
		}

		// load Link
		if((mu==0)&&(updown==1))
		{
//			COPY_GLOBALLINKTYPE( GlobalLinkType, GlobalLinkType2, 1 );
			GlobalLinkType2 globalLink(Udown, site, mu );
			localLink = globalLink;

			SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

			// write back Link
			globalLink = localLink;
		}
		else
		{
			GlobalLinkType globalLink( Uup, site, mu );
			localLink = globalLink;

			SubgroupIterator<LocalLinkType::PARAMTYPE::NC>::iterate( *this );

			// write back Link
			globalLink = localLink;
		}
	}


	__device__ inline void subgroupStep( lat_group_index_t iSub, lat_group_index_t jSub )
	{
		extern __shared__ typename LocalLinkType::PARAMTYPE::REALTYPE shA[]; // define size in kernel call (size needs to be 4*SitesPerBlock)!
		initializeSharedMemory( shA, GaugeType::SharedArraySize );

		LocalLink<SU2Vector4<typename LocalLinkType::PARAMTYPE::REALTYPE> > quaternion = localLink.getSU2Subgroup( iSub, jSub );
		typename Real4<typename LocalLinkType::PARAMTYPE::REALTYPE>::VECTORTYPE& q = quaternion[0];

		GaugeType::gatherInfo( shA, q, id, mu, updown, SitesPerBlock );

		// calc update
		__syncthreads();
		if( mu == 0 )
		{
			if( updown == 0 )
			{
				GaugeType::prepareUpdate( shA, id, SitesPerBlock );
				algorithm.calculateUpdate( shA, id, SitesPerBlock );
			}
		}
		__syncthreads();

		GaugeType::collectUpdate( shA, q, id, mu, updown, SitesPerBlock );
		__syncthreads(); // this is necessary otherwise the next subgroup already overwrites shA!
		quaternion.set( 0, q );

		// update links
		if( updown == 0 )
		{
			localLink.leftSubgroupMult( quaternion, iSub, jSub );
		}
		else
		{
			quaternion.hermitian();
			localLink.rightSubgroupMult( quaternion, iSub, jSub );
		}
	}
private:
	static const int threadsPerSite = 8;
	Algorithm algorithm;
	const int id;
	const int mu;
	bool updown;
	LocalLinkType localLink;

	__device__ inline void initializeSharedMemory( typename LocalLinkType::PARAMTYPE::REALTYPE* shA, int sharedArraySize )
	{
		if( mu == 0 )
		{
			if( updown == 0 )
			{
				for( int i = 0; i < sharedArraySize; i++ )
				{
					shA[id+i*SitesPerBlock] = 0;
				}
			}
		}
		__syncthreads();
	}
};

#endif /* GAUGEFIXINGSUBGROUPSTEP_H_ */
