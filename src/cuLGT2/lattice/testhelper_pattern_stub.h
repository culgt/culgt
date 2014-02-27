/**
 *
 *  Created on: Feb 27, 2014
 *      Author: vogt
 */

#ifndef TESTHELPER_PATTERN_STUB_H_
#define TESTHELPER_PATTERN_STUB_H_

class SiteStub
{
public:
	static const int Ndim=4;
};

template<typename T> class ParamStub
{
public:
	typedef T TYPE;
	static const int SIZE=18;
};

template<typename T> class PatternStub
{
public:
	typedef ParamStub<T> PARAMTYPE;
	typedef SiteStub SITETYPE;
};


#endif /* TESTHELPER_PATTERN_STUB_H_ */
