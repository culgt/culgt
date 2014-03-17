/**
 * LatticeDimension.h
 *
 *  Created on: Mar 6, 2014
 *      Author: vogt
 */

#ifndef LATTICEDIMENSION_H_
#define LATTICEDIMENSION_H_

#include "../common/culgt_typedefs.h"
#include "../cudacommon/cuda_host_device.h"

namespace culgt
{

template<int Ndim> class LatticeDimension
{
private:
	lat_coord_t size[Ndim];
	lat_dim_t latSize;

public:
	static const int NDIM = Ndim;

	CUDA_HOST_DEVICE inline LatticeDimension( const int size[Ndim] )
	{
		latSize = 1;
		for( int i = 0; i < Ndim; i++ )
		{
			this->size[i] = size[i];
			latSize *= size[i];
		}
	}

	/**
	 * TODO: should use variadic templates to ensure correct size constructor (but not supported in CUDA 6.0RC)
	 * @param dir0
	 * @param dir1
	 * @param dir2
	 * @param dir3
	 */
	CUDA_HOST_DEVICE inline LatticeDimension( lat_coord_t dir0 = 0, lat_coord_t dir1 = 0, lat_coord_t dir2 = 0, lat_coord_t dir3 = 0)
	{
		if( Ndim > 0 )
			size[0] = dir0;
		if( Ndim > 1 )
			size[1] = dir1;
		if( Ndim > 2 )
			size[2] = dir2;
		if( Ndim > 3 )
			size[3] = dir3;

		latSize = 1;
		for( int i = 0; i < Ndim; i++ )
		{
			latSize *= size[i];
		}
	}



	CUDA_HOST_DEVICE inline lat_coord_t getDimension( const lat_dim_t i ) const
	{
		return size[i];
	}

	CUDA_HOST_DEVICE inline lat_dim_t getSize() const
	{
		return latSize;
	}

	/**
	 * This allows the use of LatticeDimension in std::map as key.
	 */
	CUDA_HOST_DEVICE inline bool operator==( const LatticeDimension<Ndim>& rhs )
	{
		for( int i = 0; i < Ndim; i++ )
		{
			if( getDimension(i) != rhs.getDimension(i) )
				return false;
		}
		return true; // equal
	}

	CUDA_HOST_DEVICE inline bool operator!=( const LatticeDimension<Ndim>& rhs ) {return !( this == rhs);}
};


struct LatticeDimensionCompare
{
	template<typename LatticeDimensionType> bool operator()( const LatticeDimensionType& lhs, const LatticeDimensionType& rhs ) const
	{
		for( int i = 0; i < LatticeDimensionType::NDIM; i++ )
		{
			if( lhs.getDimension(i) < rhs.getDimension(i) )
				return true;
			else if( lhs.getDimension(i) > rhs.getDimension(i) )
				return false;
			// else the extent in i-direction is same, test next...
		}
		return false; // equal
	}
};

//template<int Ndim> CUDA_HOST_DEVICE inline bool operator<( const LatticeDimension<Ndim>& lhs, const LatticeDimension<Ndim>& rhs )
//{
//	for( int i = 0; i < Ndim; i++ )
//	{
//		if( lhs.getDimension(i) < rhs.getDimension(i) )
//			return true;
//		else if( lhs.getDimension(i) > rhs.getDimension(i) )
//			return false;
//		// else the extent in i-direction is same, test next...
//	}
//	return false; // equal
//}
//

}

#endif /* LATTICEDIMENSION_H_ */






