/**
 * PlaquetteAverage.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTEAVERAGE_H_
#define PLAQUETTEAVERAGE_H_

#include "../lattice/LatticeDimension.h"
#include "../lattice/LocalLink.h"
#include "../lattice/parameterization_types/SUNRealFull.h"
#include "PlaquetteAverageCudaHelper.h"
#include "../cuLGT1legacy/Reduction.hxx"

namespace culgt
{

template<typename PatternType, typename LocalLinkType = LocalLink<SUNRealFull<PatternType::PARAMTYPE::NC, typename PatternType::PARAMTYPE::TYPE> > > class PlaquetteAverage
{
public:
	typedef typename PatternType::PARAMTYPE::TYPE T;
	typedef typename PatternType::PARAMTYPE::REALTYPE REALT;

	PlaquetteAverage( T* U, LatticeDimension<PatternType::SITETYPE::Ndim> dim ) : dim(dim), U(U)
	{
		cudaMalloc( (void**)&devPtr, sizeof( REALT )*dim.getSize() );
	};

	void calculatePlaquettes()
	{
		PlaquetteAverageCudaHelper<PatternType,LocalLinkType>::calculatePlaquettes( U, devPtr, dim );
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
		cudaFree( devPtr );
	}

private:
	LatticeDimension<PatternType::SITETYPE::Ndim> dim;
	T* U;
	REALT* devPtr;
};

}

#endif /* PLAQUETTEAVERAGE_H_ */
