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

template<int Ndim> class SetSizeHelper
{
public:
	template<int dir> CUDA_HOST_DEVICE static inline void set( lat_coord_t* array, lat_coord_t size )
	{
		if( dir < Ndim )
			array[dir] = size;
	}
};


template<int Ndim> class LatticeDimension
{
private:
	lat_coord_t size[Ndim];
	lat_dim_t latSize;

	CUDA_HOST_DEVICE inline void updateLatticeSize()
	{
		latSize = 1;
		for( int i = 0; i < Ndim; i++ )
		{
			latSize *= size[i];
		}
	}
public:
	static const int NDIM = Ndim;

	CUDA_HOST_DEVICE inline LatticeDimension( const int size[Ndim] )
	{
		for( int i = 0; i < Ndim; i++ )
		{
			this->size[i] = size[i];
		}
		updateLatticeSize();
	}

	CUDA_HOST_DEVICE inline LatticeDimension( const LatticeDimension<Ndim>& src )
	{
		for( int i = 0; i < Ndim; i++ )
		{
			this->size[i] = src.getDimension(i);
		}
		updateLatticeSize();
	}

	/**
	 * TODO use variadic templates:
	 * - should use variadic templates to ensure correct size constructor (but not supported in CUDA 6.0RC)
	 * - naive code raises a warning "subscript out of range" if Ndim < 4, SetSizeHelper avoids this but i don't like it (and not well tested)
	 * @param dir0
	 * @param dir1
	 * @param dir2
	 * @param dir3
	 */
	CUDA_HOST_DEVICE inline LatticeDimension( lat_coord_t dir0 = 0, lat_coord_t dir1 = 0, lat_coord_t dir2 = 0, lat_coord_t dir3 = 0)
	{
		SetSizeHelper<Ndim>::template set<0>( size, dir0 );
		SetSizeHelper<Ndim>::template set<1>( size, dir1 );
		SetSizeHelper<Ndim>::template set<2>( size, dir2 );
		SetSizeHelper<Ndim>::template set<3>( size, dir3 );

//#pragma GCC diagnostic ignored "-Warray-bounds" // DOES NOT WORK

//		if( Ndim > 0 )
//			size[0] = dir0;
//		if( Ndim > 1 )
//			size[1] = dir1;
//		if( Ndim > 2 )
//			size[2] = dir2;
//		if( Ndim > 3 )
//			size[3] = dir3;

//#pragma GCC diagnostic pop

		updateLatticeSize();
	}



	CUDA_HOST_DEVICE inline lat_coord_t getDimension( const lat_dim_t i ) const
	{
		return size[i];
	}

	CUDA_HOST_DEVICE inline void setDimension( const lat_dim_t i, const lat_coord_t sizeInI )
	{
		size[i] = sizeInI;
		updateLatticeSize();
	}


	CUDA_HOST_DEVICE inline lat_dim_t getSize() const
	{
		return latSize;
	}

	CUDA_HOST_DEVICE inline LatticeDimension<Ndim> getDimensionTimeslice() const
	{
		LatticeDimension<Ndim> result( *this );
		result.setDimension( 0, 1 );
		return result;
	}

	CUDA_HOST_DEVICE inline LatticeDimension<Ndim-1> getReducedDimension( int dir ) const
	{
		LatticeDimension<Ndim-1> result;
		int counter = 0;
		for( int i = 0; i < Ndim; i++ )
		{
			if( i != dir )
			{
				result.setDimension( counter, this->getDimension( i ) );
				counter++;
			}
		}
		return result;
	}

	CUDA_HOST_DEVICE inline LatticeDimension<Ndim> getFFTComplexDimension() const
	{
		LatticeDimension<Ndim> result(*this);
		result.setDimension( Ndim-1, getDimension(Ndim-1)/2+1 );
		return result;
	}


	/**
	 * I don't know if this is a wise choice or if we should keep a local variable for latSizeTimeslice (but probably the compiler is smart enough to do so).
	 */
	CUDA_HOST_DEVICE inline lat_dim_t getSizeTimeslice() const
	{
		return latSize/size[0];
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

/**
 * Used in SiteNeighbourTableManager to allow LatticeDimension objects as keys in std::map.
 */
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


}

#endif /* LATTICEDIMENSION_H_ */
