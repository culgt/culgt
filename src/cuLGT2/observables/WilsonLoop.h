/**
 * Plaquette.h
 *
 *  Created on: Mar 13, 2014
 *      Author: vogt
 */

#ifndef WILSONLOOP_H_
#define WILSONLOOP_H_

#include "../lattice/LatticeDimension.h"
#include "../lattice/LocalLink.h"
#include "../lattice/GlobalLink.h"
#include "../lattice/parameterization_types/SUNRealFull.h"
#include "../cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> > > class WilsonLoop
{
public:
	typedef typename PatternType::PARAMTYPE::TYPE T;

	CUDA_HOST_DEVICE inline WilsonLoop( T* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim ) : U(U), dim(dim)
	{
	}

	CUDA_HOST_DEVICE inline LocalLinkType getWilsonLoop( typename PatternType::SITETYPE site, lat_dim_t dir1, lat_coord_t length1, lat_dim_t dir2, lat_coord_t length2 )
	{
		LocalLinkType loop;
		loop.identity();

		for( int i = 0; i < length1; i++ )
		{
			loop *= getLink( site, dir1 );
			site.setNeighbour( dir1, true );
		}

		for( int i = 0; i < length2; i++ )
		{
			loop *= getLink( site, dir2 );
			site.setNeighbour( dir2, true );
		}

		for( int i = 0; i < length1; i++ )
		{
			site.setNeighbour( dir1, false );
			loop *= getLinkHermitian( site, dir1 );
		}

		for( int i = 0; i < length2; i++ )
		{
			site.setNeighbour( dir2, false );
			loop *= getLinkHermitian( site, dir2 );
		}

		return loop;
	}

	CUDA_HOST_DEVICE inline T getReTraceWilsonLoopNormalized( typename PatternType::SITETYPE site, lat_dim_t dir1, lat_coord_t length1, lat_dim_t dir2, lat_coord_t length2 )
	{
		return getWilsonLoop( site, dir1, length1, dir2, length2 ).reTrace() / (T)PatternType::PARAMTYPE::NC;
	}


private:
	T* U;
	LatticeDimension<PatternType::SITETYPE::Ndim> dim;

	CUDA_HOST_DEVICE inline LocalLinkType getLink( const typename PatternType::SITETYPE& site, const lat_dim_t& mu )
	{
		LocalLinkType result;
		GlobalLink<PatternType> glob( U, site, mu );
		result = glob;
		return result;
	}

	CUDA_HOST_DEVICE inline LocalLinkType getLinkHermitian( const typename PatternType::SITETYPE& site, const lat_dim_t& mu )
	{
		return getLink( site, mu ).hermitian();
	}


};

}

#endif /* PLAQUETTE_H_ */

