/**
 * Plaquette.h
 *
 *  Created on: Mar 13, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTE_H_
#define PLAQUETTE_H_

#include "../lattice/LatticeDimension.h"
#include "../lattice/LocalLink.h"
#include "../lattice/GlobalLink.h"
#include "../lattice/parameterization_types/SUNRealFull.h"
#include "../cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> > > class Plaquette
{
public:
	typedef typename PatternType::PARAMTYPE::TYPE T;

	CUDA_HOST_DEVICE inline Plaquette( T* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim ) : U(U), dim(dim)
	{
	}

	CUDA_HOST_DEVICE inline LocalLinkType getStaple( typename PatternType::SITETYPE site, lat_dim_t mu, lat_dim_t nu )
	{
		LocalLinkType staple;

		staple = getLink( site.getNeighbour( mu, true ), nu );
		staple *= getLinkHermitian( site.getNeighbour( nu, true), mu );
		staple *= getLinkHermitian( site, nu );

		return staple;
	}

	CUDA_HOST_DEVICE inline LocalLinkType getPlaquette( typename PatternType::SITETYPE site, lat_dim_t mu, lat_dim_t nu )
	{
		LocalLinkType temp = getLink( site, mu );
		temp *= getStaple( site, mu, nu );
		return temp;
	};

	CUDA_HOST_DEVICE inline T getReTracePlaquetteNormalized( typename PatternType::SITETYPE site, lat_dim_t mu, lat_dim_t nu )
	{
		return getPlaquette( site, mu, nu ).reTrace() / (T)PatternType::PARAMTYPE::NC;
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

