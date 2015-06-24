/**
 *
 * Lattice Site using a Nd-dimensional array.
 *
 * This implementation uses a Nd-dimensional array to store the lattice coordinates of the current site.
 * We should use the "SiteIndex" implementation whenever we go for minimal register usage.
 *
 * Possible Optimizations:
 *  - We use on the fly neighbour calculation here. Check if a table is a more favorable choice.
 *  - See "SiteIndex".
 */

#ifndef SITECOORD_H_
#define SITECOORD_H_

#include "ParityType.h"
#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"
#include "lattice/LatticeDimension.h"
#include <assert.h>

namespace culgt
{

/**
 * The template parameter "ParityType par" defines normal indexing (par==NO_SPLIT), full parity splitting (par==FULL_SPLIT) or splitting
 * only within a timeslice (par == TIMESLICE_SPLIT)
 * Parity split indexing means that the first half of the indices are the even sites ((x+y+z+...)%2 == 0) and the second half the odd sites.
 * It has to be used for CUDA implementations to ensure coalesced memory reads for neighbouring (concerning its index) array elements.
 *
 * @param template Nd: Dimension of the lattice
 * @param template par
 */
template<lat_dim_t Nd, ParityType par> class SiteCoord
{
public:
	lat_coord_t size[Nd];
	static const lat_dim_t NDIM = Nd;
	static const ParityType PARITYTYPE = par;
	lat_coord_t site[Nd];

	CUDA_HOST_DEVICE inline SiteCoord( const lat_coord_t size[Nd], lat_index_t* nn=NULL )
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = size[i];
		}
	}

	CUDA_HOST_DEVICE inline SiteCoord( const culgt::LatticeDimension<Nd>& dim, lat_index_t* nn=NULL )
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = dim.getDimension( i );
		}
	}

	CUDA_HOST_DEVICE inline SiteCoord( const culgt::LatticeDimension<Nd> dim, NeigbourTableType type )
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = dim.getDimension( i );
		}
	}

	CUDA_HOST_DEVICE inline SiteCoord( const SiteCoord<Nd,par> &s )
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = s.size[i];
			this->site[i] = s.site[i];
		}
	}

	CUDA_HOST_DEVICE inline void calculateNeighbourTable( lat_index_t* nn ) {}; // does not use neighbourtable (for compatibility with SiteIndex);

	CUDA_HOST_DEVICE inline lat_coord_t& operator[](lat_dim_t i)
	{
		return site[i];
	}

	CUDA_HOST_DEVICE inline lat_coord_t getCoord( const lat_dim_t i) const
	{
		return site[i];
	}

	/**
	 * Returns the size of the lattice in direction i
	 * @return size of timeslice
	 */
	CUDA_HOST_DEVICE inline lat_coord_t getLatticeSizeDirection( lat_dim_t i )
	{
		return size[i];
	}

	CUDA_HOST_DEVICE inline lat_index_t getIndex() const
	{
		if( par == FULL_SPLIT )
		{
			lat_index_t parity = 0;
			for(  lat_dim_t i = 0; i < Nd; i++ )
			{
				parity += site[i];
			}

			lat_index_t index = 0;
			for( lat_dim_t i = 0; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}

			if( parity % 2 == 0 )
			{
				return index / 2;
			}
			else
			{
				return index / 2 + getSize()/2;
			}
		}
		else if( par == TIMESLICE_SPLIT )
		{
			// TODO check this for correctness
			lat_index_t parity = 0;
			for(  lat_dim_t i = 1; i < Nd; i++ )
			{
				parity += site[i];
			}

			lat_index_t index = 0;
			for( lat_dim_t i = 1; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}

			if( parity % 2 == 0 )
			{
				index /= 2;
			}
			else
			{
				index = index / 2 + getSizeTimeslice()/2;
			}
			return index + site[0] * getSizeTimeslice();
		}
		else // NO_SPLIT
		{
			lat_index_t index = 0;
			for( lat_dim_t i = 0; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}
			return index;
		}
	}

	CUDA_HOST_DEVICE inline lat_index_t getIndexNonSplit() const
	{
		lat_index_t index = 0;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			index *= size[i];
			index += site[i];
		}
		return index;
	}

	CUDA_HOST_DEVICE inline void setIndex( lat_index_t latticeIndex )
	{
		// TODO simplify this
		bool curParity;
		if( par == FULL_SPLIT )
		{
			if( latticeIndex >= getSize() / 2 ) // odd parity
			{
				latticeIndex -= getSize() / 2;
				latticeIndex *= 2;
				curParity = 1;
			}
			else // even parity
			{
				latticeIndex *= 2;
				curParity = 0;
			}

			lat_index_t parity = 0;
			for( lat_dim_t i = Nd-1; i >= 0; i-- )
			{
				site[i] = latticeIndex % size[i];
				parity += site[i];
				latticeIndex /= size[i];
			}
			if( parity % 2 != curParity )
				site[Nd-1]++;
		}
		if( par == NO_SPLIT )
		{
			lat_index_t parity = 0;
			for( lat_dim_t i = Nd-1; i >= 0; i-- )
			{
				site[i] = latticeIndex % size[i];
				parity += site[i];
				latticeIndex /= size[i];
			}

		}

		// TODO test this!!!
		if( par == TIMESLICE_SPLIT )
		{
			site[0] = latticeIndex / getSizeTimeslice();
			lat_index_t tempIndex = latticeIndex % getSizeTimeslice(); // remove temporal index
			if( tempIndex >= getSizeTimeslice() / 2 ) // odd parity
			{
				tempIndex -= getSizeTimeslice() / 2;
				tempIndex *= 2;
				curParity = 1;
			}
			else // even parity
			{
				tempIndex *= 2;
				curParity = 0;
			}

			lat_index_t parity = 0;
			for( lat_dim_t i = Nd-1; i >= 1; i-- )
			{
				site[i] = tempIndex % size[i];
				parity += site[i];
				tempIndex /= size[i];
			}
			if(parity % 2 != curParity )
				site[Nd-1]++;

		}


	}


	CUDA_HOST_DEVICE inline void setIndexFromParitySplitOrder( lat_index_t latticeIndex )
	{
		// TODO BUG! Compare to setLatticeIndex()!!!

		if( latticeIndex >= getSize() / 2 )
		{
			latticeIndex -= getSize() / 2;
			latticeIndex *= 2;
			latticeIndex++;
		}
		else
		{
			latticeIndex *= 2;
		}

		for( lat_dim_t i = Nd-1; i >= 0; i-- )
		{
			site[i] = latticeIndex % size[i];
			latticeIndex /= size[i];
		}
	}

	CUDA_HOST_DEVICE inline void setIndexFromNonParitySplitOrder( lat_index_t latticeIndex )
	{
		for( lat_dim_t i = Nd-1; i >= 0; i-- )
		{
			site[i] = latticeIndex % size[i];
			latticeIndex /= size[i];
		}
	}


	CUDA_HOST_DEVICE inline lat_index_t getSize() const
	{
		int tmp = 1;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			tmp *= size[i];
		}
		return tmp;
	}

	/**
	 * Index within a timeslice
	 */
	CUDA_HOST_DEVICE inline lat_index_t getIndexTimeslice() const
	{
		if( par == FULL_SPLIT )
		{
			// TODO this looks wrong for FULL_SPLIT. Looks more like what TIMESLICE_SPLIT wants to do (therefore I copy it)...
			// maybe this is good for backward compatibility.
			// as I see it atm CoulombGauge should use TIMESLICE_SPLIT but uses this
			// OR BETTER: REMOVE THIS and check that all apps use the intuitive TIMESLICE_SPLIT!!!
//					assert(false);

			lat_index_t parity = 0;
			for(  lat_dim_t i = 1; i < Nd; i++ )
			{
				parity += site[i];
			}

			lat_index_t index = 0;
			for( lat_dim_t i = 1; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}

			if( parity % 2 == 0 )
			{
				return index / 2;
			}
			else
			{
				return index / 2 + getSizeTimeslice()/2;
			}
		}
		else if( par == TIMESLICE_SPLIT )
		{
			lat_index_t parity = 0;
			for(  lat_dim_t i = 1; i < Nd; i++ )
			{
				parity += site[i];
			}

			lat_index_t index = 0;
			for( lat_dim_t i = 1; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}

			if( parity % 2 == 0 )
			{
				return index / 2;
			}
			else
			{
				return index / 2 + getSizeTimeslice()/2;
			}
		}
		else // NO_SPLIT
		{
			lat_index_t index = 0;
			for( lat_dim_t i = 1; i < Nd; i++ )
			{
				index *= size[i];
				index += site[i];
			}
			return index;
		}
	}

	CUDA_HOST_DEVICE inline lat_index_t getSizeTimeslice() const
	{
		int tmp = 1;
		for( lat_dim_t i = 1; i < Nd; i++ )
		{
			tmp *= size[i];
		}
		return tmp;
	}


	/**
	 * If you think this for-loop looks strange: This is a way the compiler can already index the array! -> The site[Ndim] array is not placed in (CUDA) memory
	 * TODO maybe check if newer compiler can do it directly
	 */
	CUDA_HOST_DEVICE inline void setNeighbour( lat_dim_t direction, bool up )
	{
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			if( direction == i )
			{
				if( up )
					site[i]++;
				else
					site[i]--;

				if( site[i] >= size[i] )
				{
					site[i] -= size[i];
				}
				else
				{
					if( site[i] < 0 )
					{
						site[i] += size[i];
					}
				}
			}
		}
	}

	/**
	 * calculates neighbour (t,x,y,z) + (direction[0],direction[1],...), i.e. moves with vector direction on the lattice
	 *
	 * If you think this for-loop looks strange: This is a way the compiler can already index the array! -> The site[Ndim] array is not placed in (CUDA) memory
	 */
	CUDA_HOST_DEVICE inline void setNeighbour( lat_coord_t* direction )
	{
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			{
				site[i] += direction[i];

				if( site[i] >= size[i] )
				{
					site[i] -= size[i];
				}
				else
				{
					if( site[i] < 0 )
					{
						site[i] += size[i];
					}
				}
			}
		}
	}

	CUDA_HOST_DEVICE inline SiteCoord<Nd, par> getNeighbour( lat_dim_t direction, bool up )
	{
		SiteCoord<Nd, par> site(*this);
		site.setNeighbour( direction, up );
		return site;
	}

};

}

#endif /* SITECOORD_HXX_ */
