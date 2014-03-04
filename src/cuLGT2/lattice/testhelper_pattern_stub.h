/**
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef TESTHELPER_PATTERN_STUB_H_
#define TESTHELPER_PATTERN_STUB_H_

template<int NDIM = 4> class SiteStub
{
public:
	static const int Ndim=NDIM;
};

template<typename T, int TNC=3> class ParamStub
{
public:
	typedef T TYPE;
	static const int NC = TNC;
	static const int SIZE=18;
};

template<typename T, int SiteStubNdim=4, int ParamStubNc=3> class PatternStub
{
public:
	typedef ParamStub<T, ParamStubNc> PARAMTYPE;
	typedef SiteStub<SiteStubNdim> SITETYPE;
};


#endif /* TESTHELPER_PATTERN_STUB_H_ */
