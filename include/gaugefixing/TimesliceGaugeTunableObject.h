/**
 *
 *  Created on: Apr 30, 2014
 *      Author: vogt
 */

#ifndef TIMESLICEGAUGETUNABLEOBJECT_H_
#define TIMESLICEGAUGETUNABLEOBJECT_H_

#include "lattice/GaugeConfigurationCudaHelper.h"
#include "util/performance/TunableObject.h"

namespace culgt
{

template<typename GlobalLinkType, typename LocalLinkType> class TimesliceGaugeTunableObject: public TunableObject
{
public:
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::TYPE T;
	typedef typename GlobalLinkType::PATTERNTYPE::PARAMTYPE::REALTYPE REALT;
	TimesliceGaugeTunableObject( T** Ut, T** UtDown, LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dimTimeslice, long seed ) : Ut(Ut), UtDown(UtDown), dimTimeslice(dimTimeslice), seed(seed)
	{
	}

	void preTuneAction()
	{
		GaugeConfigurationCudaHelper<T>::template setHot<typename GlobalLinkType::PATTERNTYPE,PhiloxWrapper<REALT> >( *Ut, dimTimeslice, seed, PhiloxWrapper<REALT>::getNextCounter() );
		GaugeConfigurationCudaHelper<T>::template setHot<typename GlobalLinkType::PATTERNTYPE,PhiloxWrapper<REALT> >( *UtDown, dimTimeslice, seed, PhiloxWrapper<REALT>::getNextCounter() );
	}

	void setSeed( long seed )
	{
		this->seed = seed;
	}

protected:
	T** Ut;
	T** UtDown;
	LatticeDimension<GlobalLinkType::PATTERNTYPE::SITETYPE::NDIM> dimTimeslice;
	long seed;
};




}

#endif
