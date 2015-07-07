/**
 * PlaquetteAverage.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTEAVERAGE_H_
#define PLAQUETTEAVERAGE_H_

#include "cudacommon/DeviceMemoryManager.h"
#include "lattice/LatticeDimension.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "PlaquetteAverageCudaHelper.h"
#include "util/reduction/Reduction.h"

namespace culgt
{

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> >, PlaquetteType MyPlaquetteType = PLAQUETTE_FULL > class PlaquetteAverage
{
public:
	typedef typename PatternType::PARAMTYPE::TYPE T;
	typedef typename PatternType::PARAMTYPE::REALTYPE REALT;

	PlaquetteAverage( T* U, LatticeDimension<PatternType::SITETYPE::NDIM> dim ) : dim(dim), U(U)
	{
		DeviceMemoryManager::malloc( &devPtr, sizeof( REALT )*dim.getSize(), "local plaquette" );
	}

	void calculatePlaquettes()
	{
		PlaquetteAverageCudaHelper<PatternType,LocalLinkType,MyPlaquetteType>::calculatePlaquettes( U, devPtr, dim );
	}

	REALT getPlaquette()
	{
		calculatePlaquettes();
		Reduction<REALT> reducer(dim.getSize());

		return reducer.reduceAll( devPtr )/(REALT)dim.getSize();
	}

	REALT* getDevicePointer() const
	{
		return devPtr;
	}

	~PlaquetteAverage()
	{
		DeviceMemoryManager::free( devPtr );
	}

private:
	LatticeDimension<PatternType::SITETYPE::NDIM> dim;
	T* U;
	REALT* devPtr;
};

}

#endif /* PLAQUETTEAVERAGE_H_ */
