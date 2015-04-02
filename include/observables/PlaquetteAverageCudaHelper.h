/**
 * PlaquetteAverage.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTEAVERAGECUDAHELPER_H_
#define PLAQUETTEAVERAGECUDAHELPER_H_

#include "lattice/LatticeDimension.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/KernelSetup.h"
#include "cudacommon/cuda_error.h"
#include "Plaquette.h"
#include "lattice/site_indexing/SiteNeighbourTableManager.h"

namespace culgt
{

enum PlaquetteType { PLAQUETTE_FULL, PLAQUETTE_SPATIAL, PLAQUETTE_TEMPORAL };

namespace PlaquetteAverageCudaHelperKernel
{
	template<typename PatternType, typename LocalLinkType, PlaquetteType MyPlaquetteType> __global__ void kernelCalculatePlaquettes( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* plaquetteArray,  LatticeDimension<PatternType::SITETYPE::NDIM> dim, lat_index_t* nn )
	{
		typedef typename PatternType::PARAMTYPE::REALTYPE T;
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site(dim, nn );
		site.setIndex( index );

		Plaquette<PatternType,LocalLinkType> p( U, dim );

		T temp = 0;

		int startmu = 0;
		int endmu = PatternType::SITETYPE::NDIM;

		if( MyPlaquetteType == PLAQUETTE_SPATIAL ) startmu = 1;
		if( MyPlaquetteType == PLAQUETTE_TEMPORAL ) endmu = 1;

		for( int mu = startmu; mu < endmu; mu++ )
			for( int nu = mu+1; nu < PatternType::SITETYPE::NDIM; nu++ )
			{
				temp += p.getReTracePlaquetteNormalized( site, mu, nu );

			}

		lat_coord_t norm;
		if( MyPlaquetteType == PLAQUETTE_SPATIAL )
		{
			norm = (PatternType::SITETYPE::NDIM-1)*(PatternType::SITETYPE::NDIM-2)/2;
		}
		else if( MyPlaquetteType == PLAQUETTE_TEMPORAL )
		{
			norm = (PatternType::SITETYPE::NDIM-1);
		}
		else
		{
			norm = PatternType::SITETYPE::NDIM*(PatternType::SITETYPE::NDIM-1)/2;
		}
		plaquetteArray[index] = temp/(T)norm;

	};
}

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> >, PlaquetteType MyPlaquetteType = PLAQUETTE_FULL> class PlaquetteAverageCudaHelper
{
public:

	static void calculatePlaquettes( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* plaquetteArray,  LatticeDimension<PatternType::SITETYPE::NDIM> dim )
	{
		KernelSetup<PatternType::SITETYPE::NDIM> setup( dim );
		PlaquetteAverageCudaHelperKernel::kernelCalculatePlaquettes<PatternType,LocalLinkType,MyPlaquetteType><<<setup.getGridSize(),setup.getBlockSize()>>>( U, plaquetteArray, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
		CUDA_LAST_ERROR( "kernelCalculatePlaquettes" );
	}

};

}

#endif
