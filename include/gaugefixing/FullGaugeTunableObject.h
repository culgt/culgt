/**
 *  Created on: Apr 30, 2014
 *      Author: vogt
 */

#ifndef FULLGAUGETUNABLEOBJECT_H_
#define FULLGAUGETUNABLEOBJECT_H_

#include "../util/performance/TunableObject.h"
#include "../util/rng/PhiloxWrapper.h"

namespace culgt
{

template<typename GlobalLinkType, typename LocalLinkType> class FullGaugeTunableObject: public TunableObject
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	FullGaugeTunableObject( T** U, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dim, long seed ) : U(U), dim(dim), seed(seed)
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
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dim;
	long seed;
};




}

#endif
