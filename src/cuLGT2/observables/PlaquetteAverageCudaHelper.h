/**
 * PlaquetteAverage.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTEAVERAGECUDAHELPER_H_
#define PLAQUETTEAVERAGECUDAHELPER_H_

#include "../lattice/LatticeDimension.h"
#include "../lattice/LocalLink.h"
#include "../lattice/parameterization_types/SUNRealFull.h"
#include "../cudacommon/cuda_error.h"
#include "../lattice/KernelSetup.h"
#include "Plaquette.h"
#include "../cuLGT1legacy/SiteNeighbourTableManager.h"

namespace culgt
{

namespace PlaquetteAverageCudaHelperKernel
{
	template<typename PatternType, typename LocalLinkType> __global__ void kernelCalculatePlaquettes( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* plaquetteArray,  LatticeDimension<PatternType::SITETYPE::Ndim> dim, lat_index_t* nn )
	{
		typedef typename PatternType::PARAMTYPE::REALTYPE T;
		int index = blockIdx.x * blockDim.x + threadIdx.x;

		VERIFY_LATTICE_SIZE( dim, index );

		typename PatternType::SITETYPE site(dim, nn );
		site.setLatticeIndex( index );

		Plaquette<PatternType,LocalLinkType> p( U, dim );

		T temp = 0;
		for( int mu = 0; mu < PatternType::SITETYPE::Ndim; mu++ )
			for( int nu = mu+1; nu < PatternType::SITETYPE::Ndim; nu++ )
			{
				temp += p.getReTracePlaquetteNormalized( site, mu, nu );

			}
		lat_coord_t norm = PatternType::SITETYPE::Ndim*(PatternType::SITETYPE::Ndim-1)/2;
		plaquetteArray[index] = temp/(T)norm;

	};
}

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> > > class PlaquetteAverageCudaHelper
{
public:

	static void calculatePlaquettes( typename PatternType::PARAMTYPE::TYPE* U, typename PatternType::PARAMTYPE::REALTYPE* plaquetteArray,  LatticeDimension<PatternType::SITETYPE::Ndim> dim )
	{
		KernelSetup<PatternType::SITETYPE::Ndim> setup( dim );
		PlaquetteAverageCudaHelperKernel::kernelCalculatePlaquettes<PatternType,LocalLinkType><<<setup.getGridSize(),setup.getBlockSize()>>>( U, plaquetteArray, dim, SiteNeighbourTableManager<typename PatternType::SITETYPE>::getDevicePointer( dim ) );
		CUDA_LAST_ERROR( "kernelCalculatePlaquettes" );
	}

};

}

#endif
