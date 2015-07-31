/**
 * test_PatternMocks.h
 *
 *  Created on: Feb 20, 2014
 *      Author: vogt
 */

#ifndef TEST_PATTERNMOCKS_H_
#define TEST_PATTERNMOCKS_H_
#include "lattice/site_indexing/ParityType.h"
#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"

using namespace culgt;

template<int TNdim=4, ParityType Parity=NO_SPLIT> class SiteTypeMock
{
public:
	static const ParityType PARITYTYPE = Parity;
	static const int NDIM = TNdim;
	const int size;
	const int sizeTimeslice;
	lat_index_t index;
	lat_index_t indexTimeslice;
	CUDA_HOST_DEVICE SiteTypeMock( int size, int sizeTimeslice, int index, int indexTimeslice ) : size(size), sizeTimeslice(sizeTimeslice), index(index), indexTimeslice(indexTimeslice) {};
	CUDA_HOST_DEVICE SiteTypeMock( int size, int index ) : size(size), sizeTimeslice(0), index(index), indexTimeslice(0) {};
	CUDA_HOST_DEVICE SiteTypeMock( int index ) : size(0), sizeTimeslice(0), index(index), indexTimeslice(0) {};
	CUDA_HOST_DEVICE SiteTypeMock() : size(0), sizeTimeslice(0), index(0), indexTimeslice(0) {};
	CUDA_HOST_DEVICE void setIndex( int index )
	{
		this->index = index;
	}
	CUDA_HOST_DEVICE lat_index_t getIndex() const
 	{
		return index;
	}
	CUDA_HOST_DEVICE lat_index_t getIndexTimeslice() const
	{
		return indexTimeslice;
	}
	CUDA_HOST_DEVICE lat_index_t getSize() const
	{
		return size;
	}
	CUDA_HOST_DEVICE lat_index_t getSizeTimeslice() const
	{
		return sizeTimeslice;
	}

	/**
	 * always return index of timeslice
	 * @param mu
	 * @return
	 */
	CUDA_HOST_DEVICE int getCoord( const int mu) const
	{
		return index/sizeTimeslice;
	}
	CUDA_HOST_DEVICE int operator[]( const int mu ) const
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
