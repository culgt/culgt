/**
 *  Created on: Apr 30, 2014
 *      Author: vogt
 */

#ifndef LANDAUGAUGETUNABLEOBJECT_H_
#define LANDAUGAUGETUNABLEOBJECT_H_

#include "../util/performance/TunableObject.h"
#include "../util/rng/PhiloxWrapper.h"

namespace culgt
{

template<typename GlobalLinkType, typename LocalLinkType> class LandauGaugeTunableObject: public TunableObject
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	LandauGaugeTunableObject( T** U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim, long seed ) : U(U), dim(dim), seed(seed)
	{
	}

	void preTuneAction()
	{
		GaugeConfigurationCudaHelper<T>::template setHot<typename GlobalLinkType::PATTERNTYPE,PhiloxWrapper<REALT> >( *U, dim, seed, PhiloxWrapper<REALT>::getNextCounter() );
	}

	void setSeed( long seed )
	{
		this->seed = seed;
	}

protected:
	T** U;
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::Ndim> dim;
	long seed;
};




}

#endif
