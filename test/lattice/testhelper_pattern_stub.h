/**
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef TESTHELPER_PATTERN_STUB_H_
#define TESTHELPER_PATTERN_STUB_H_
#include "common/culgt_typedefs.h"
#include "lattice/LatticeDimension.h"
#include "lattice/ParameterizationMediator.h"
#include "lattice/GlobalLink.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNRealFull.h"

template<int NDIM = 4, int INDEX=0> class SiteStub
{
public:
	SiteStub(){};
	SiteStub( culgt::LatticeDimension<NDIM> dim, lat_index_t* nn = NULL ){};
	static const int Ndim=NDIM;
	void setLatticeIndex( const int i ){};
	void setLatticeIndexFromNonParitySplitOrder( const int i ){};
	int getIndex()  const { return INDEX; };
};

template<typename T, int TNC=3> class ParamStub
{
public:
	typedef T TYPE;
	typedef T REALTYPE;
	static const int NC = TNC;
	static const int SIZE=18;
};

template<typename T, int SiteStubNdim=4, int ParamStubNc=3> class PatternStub
{
public:
	typedef ParamStub<T, ParamStubNc> PARAMTYPE;
	typedef SiteStub<SiteStubNdim> SITETYPE;
	static lat_index_t getIndex( SITETYPE s, lat_dim_t mu, lat_group_index_t i ){ return 0; };
};


namespace culgt{

template<typename T, int NC, int NDIM> class ParameterizationMediator< ParamStub<T, NC>,SUNRealFull<NC, T>,GlobalLink<PatternStub<T, NDIM, NC> >,LocalLink<SUNRealFull<NC, T> > >
{
public:
	CUDA_HOST_DEVICE static void assign( GlobalLink<PatternStub<T, NDIM, NC> >& l1, const LocalLink<SUNRealFull<NC, T> >& l2 )
	{

	};
};
template<typename T, int NC, int NDIM> class ParameterizationMediator<SUNRealFull<NC, T>, ParamStub<T, NC>,LocalLink<SUNRealFull<NC, T> >,GlobalLink<PatternStub<T, NDIM, NC> > >
{
public:
	CUDA_HOST_DEVICE static void assign( LocalLink<SUNRealFull<NC, T> >& l1, const GlobalLink<PatternStub<T, NDIM, NC> >& l2 )
	{

	};
};
}

#endif /* TESTHELPER_PATTERN_STUB_H_ */
