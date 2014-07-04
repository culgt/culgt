/**
 * test_PatternMocks.h
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */

#ifndef TEST_PATTERNMOCKS_H_
#define TEST_PATTERNMOCKS_H_
#include "../../cuLGT1legacy/ParityType.h"

template<int TNdim=4, ParityType Parity=NO_SPLIT> class SiteTypeMock
{
public:
	static const ParityType PARITYTYPE = Parity;
	static const int Ndim = TNdim;
	const int size;
	const int sizeTimeslice;
	int index;
	int indexTimeslice;
	SiteTypeMock( int size, int sizeTimeslice, int index, int indexTimeslice ) : size(size), sizeTimeslice(sizeTimeslice), index(index), indexTimeslice(indexTimeslice) {};
	SiteTypeMock( int size, int index ) : size(size), sizeTimeslice(0), index(index), indexTimeslice(0) {};
	SiteTypeMock( int index ) : size(0), sizeTimeslice(0), index(index), indexTimeslice(0) {};
	void setIndex( int index )
	{
		this->index = index;
	}
	int getIndex() const
 	{
		return index;
	}
	int getIndexTimeslice() const
	{
		return indexTimeslice;
	}
	int getSize() const
	{
		return size;
	}
	int getSizeTimeslice() const
	{
		return sizeTimeslice;
	}
	/**
	 * always return index of timeslice
	 * @param mu
	 * @return
	 */
	int operator[]( const int mu ) const
	{
		return index/sizeTimeslice;
	}
};


class SiteTypeMock2
{
public:
	int index;
	int getIndex() const
 	{
		return index;
	}
	int getSize() const
	{
		return 2;
	}
};



template<int mySize, int nc=2> class ParamTypeMock
{
public:
	static const int SIZE = mySize;
	static const int NC = nc;
};

class ParamTypeMock2
{
public:
	static const int SIZE = 2;
};




#endif /* TEST_PATTERNMOCKS_H_ */
