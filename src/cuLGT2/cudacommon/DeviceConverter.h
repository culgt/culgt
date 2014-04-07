/**
 * DeviceConverter.h
 *
 *  Created on: Apr 7, 2014
 *      Author: vogt
 */

#ifndef DEVICECONVERTER_H_
#define DEVICECONVERTER_H_
#include "../lattice/KernelSetup.h"

namespace culgt
{

namespace DeviceConverterKernel
{
/**
 * TODO: if SiteType is the same we could shortcut getIndexNonSplit() to getIndex() which should be much faster for SPLIT layouts.
 * @param destU
 * @param srcU
 * @param dim
 */
	template<typename DestPattern, typename SrcPattern>__global__ void convert( DestPattern::PARAMTYPE::TYPE* destU, SrcPattern::PARAMTYPE::TYPE* srcU,LatticeDimension<DestPattern::SITETYPE::Ndim> dim )
	{

		for( int mu = 0; mu < DestPattern::SITETYPE::Ndim; mu++ )
		{
			int index = blockIdx.x * blockDim.x + threadIdx.x;

			VERIFY_LATTICE_SIZE( dim, index );

			typename SrcPattern::SITETYPE siteSrc( dim, DO_NOT_USE_NEIGHBOURS );
			siteSrc.setLatticeIndex( index );

			LocalLink<typename DestPattern::PARAMTYPE> local;
			GlobalLink<SrcPattern> src( srcU, siteSrc, mu );
			local = src;

			typename DestPattern::SITETYPE siteDest( dim, DO_NOT_USE_NEIGHBOURS );
			siteDest.setLatticeIndexFromNonParitySplitOrder( siteSrc.getIndexNonSplit() );
			GlobalLink<DestPattern> dest( destU, siteDest, mu );
			dest = local;
		}
	}
}


template<typename DestPattern, typename SrcPattern> class DeviceConverter
{
	static void convert( DestPattern::PARAMTYPE::TYPE* destU, SrcPattern::PARAMTYPE::TYPE* srcU, LatticeDimension<DestPattern::SITETYPE::Ndim> dim )
	{
		KernelSetup<DestPattern::SITETYPE::Ndim> setup( dim );

		DeviceConverterKernel::convert<DestPattern,SrcPattern><<<setup.getGridSize(),setup.getBlockSize()>>>( destU, srcU, dim );
	}


};

}

#endif /* DEVICECONVERTER_H_ */
