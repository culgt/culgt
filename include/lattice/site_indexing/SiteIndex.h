/**
 *
 * Lattice Site using a 1 dimensional index.
 *
 * This implementation of "Site" uses the 1 dimensional index to store the current site. It is favorable in CUDA implementations
 * since it uses fewer variables (exactly one) to store the site compared to "SiteCoord" which uses Ndim variables.
 * This results in fewer registers which is very important for CUDA optimization.
 * Some functions are more costly in this implementation, like operator[] (access to coordinates), other functions are much faster,
 * like getLatticeIndex() (which simply return the member variable).
 */

#ifndef SITEINDEX_H_
#define SITEINDEX_H_

#include "cudacommon/cuda_host_device.h"
#include "common/culgt_typedefs.h"
#include "lattice/LatticeDimension.h"
#include "ParityType.h"
#include <assert.h>

/**
 * See SiteCoord.hxx for a note on the parity splitting.
 *
 * @param template Nd: Dimension of the lattice
 * @param template par
 * @see SiteCoord.hxx
 */

namespace culgt
{

template<lat_dim_t Nd, ParityType par> class SiteIndex
{
private:
	lat_index_t index;
	lat_index_t latticeSize;
public:
	static const lat_dim_t NDIM = Nd;
	static const ParityType PARITYTYPE = par;
	lat_index_t* nn;
	lat_coord_t size[Nd];


	CUDA_HOST_DEVICE inline SiteIndex( const lat_coord_t size[Nd], lat_index_t* nn ) :nn(nn)
	{
		latticeSize = 1;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			this->size[i] = size[i];
			latticeSize *= size[i];
		}
		index = 0;
	}

	CUDA_HOST_DEVICE inline SiteIndex( const lat_coord_t size[Nd], NeigbourTableType type ) : nn(NULL)
	{
		latticeSize = 1;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			this->size[i] = size[i];
			latticeSize *= size[i];
		}
		index = 0;
	}

	CUDA_HOST_DEVICE inline SiteIndex( const culgt::LatticeDimension<Nd> dim, lat_index_t* nn ):nn(nn)
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = dim.getDimension( i );
		}
		latticeSize = dim.getSize();
		index = 0;
	}

	CUDA_HOST_DEVICE inline SiteIndex( const culgt::LatticeDimension<Nd> dim, NeigbourTableType type ) : nn(NULL)
	{
		for( int i = 0; i < Nd; i++ )
		{
			this->size[i] = dim.getDimension( i );
		}
		latticeSize = dim.getSize();
		index = 0;
	}

	CUDA_HOST_DEVICE inline SiteIndex( const SiteIndex<Nd,par> &s )
	{
		this->index = s.index;
	//	latticeSize = 1;
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			this->size[i] = s.size[i];
	//		latticeSize *= size[i];
		}
		this->latticeSize = s.getSize();
		this->nn = s.nn;
	}

	CUDA_HOST_DEVICE inline SiteIndex( const lat_index_t latticeSize, lat_index_t* nn ) :latticeSize(latticeSize), nn(nn), index(0)
	{
		for( lat_dim_t i = 0; i < Nd; i++ )
		{
			this->size[i] = 0;
		}
	}

	CUDA_HOST_DEVICE inline ~SiteIndex()
	{
	}

	/**
	 * TODO
	 * Returns the i-coordinate of the current site.
	 */
	CUDA_HOST_DEVICE inline lat_coord_t operator[](const lat_dim_t i) const
	{
		lat_index_t temp = index;

		if( par == FULL_SPLIT )
		{
			assert(false); // TODO
			return -1;
		}
		else if ( par == TIMESLICE_SPLIT )
		{
			if( i == 0 )
			{
				return temp / (latticeSize/size[0]);
			}
			else
			{
				assert(false); // TODO exercise for the reader
				return -1;
			}
		}
		else // NO_SPLIT
		{
			for( lat_dim_t j = Nd-1; j >= 0; j-- )
			{
				if( i == j ) return temp % size[j];
				else temp /= size[j];
			}
		}
		assert(false); // if we are here, something is wrong...
		return -1;
	}

	CUDA_HOST_DEVICE inline lat_coord_t getCoord(const lat_dim_t i) const
	{
		lat_index_t temp = index;

		if( par == FULL_SPLIT )
		{
			assert(false); // TODO
			return -1;
		}
		else if ( par == TIMESLICE_SPLIT )
		{
			if( i == 0 )
			{
				return temp / (latticeSize/size[0]);
			}
			else
			{
				assert(false); // TODO exercise for the reader
				return -1;
			}
		}
		else // NO_SPLIT
		{
			for( lat_dim_t j = Nd-1; j >= 0; j-- )
			{
				if( i == j ) return temp % size[j];
				else temp /= size[j];
			}
		}
		assert(false); // if we are here, something is wrong...
		return -1;
	}

	/**
	 * Returns the 1 dimensional index of the current site
	 * @return lattice index
	 */
	CUDA_HOST_DEVICE inline lat_index_t getIndex() const
	{
		return index;
	}

	/**
	 * Don't hate me for this
	 * @return
	 */
	CUDA_HOST_DEVICE inline lat_index_t getIndexNonSplit() const
	{
		if( par == FULL_SPLIT )
		{
			lat_index_t nonSplitIndex = index;
			bool parity = 0;
			if( nonSplitIndex >= getSize()/2 )
			{
				parity=1;
				nonSplitIndex -= getSize()/2;
			}
			nonSplitIndex *= 2;

			lat_index_t checkParity = nonSplitIndex;

			lat_coord_t curParity;
			for( lat_dim_t i = Nd-1; i >= 0; i-- )
			{
				curParity += checkParity % size[i];
				checkParity /= size[i];
			}

			curParity %= 2;
			if( !curParity && parity )
				nonSplitIndex += 1;
			else if( curParity && !parity )
				nonSplitIndex += 1;

			return nonSplitIndex;
		}
		else if ( par == TIMESLICE_SPLIT )
		{
			assert(false);
			return -1;
		}
		else
		{
			return index;
		}
	}

	/**
	 * Sets the lattice index.
	 * @param lattice index
	 * @return void
	 */
	CUDA_HOST_DEVICE inline void setIndex( lat_index_t latticeIndex )
	{
		index = latticeIndex;
	}

	/**
	 * Sets the lattice index.
	 * @param lattice index
	 * @return void
	 */
	CUDA_HOST_DEVICE inline void setIndexTimeslice( lat_index_t latticeIndex, lat_coord_t t )
	{
		if( par == FULL_SPLIT )
		{
			assert(false); // TODO
		}
		else if ( par == TIMESLICE_SPLIT )
		{
			assert(false); // TODO
		}
		else // NO_SPLIT
		{
			index = t*getSizeTimeslice() + latticeIndex;
		}
	}

	/**
	 * // TODO
	 * If par == true this simply sets the given index, if par == false the index has to be converted to the non-split index.
	 */
	CUDA_HOST_DEVICE inline void setIndexFromParitySplitOrder( lat_index_t latticeIndex )
	{
		// TODO
		assert(false);
	}

	/**
	 * // TODO
	 * If par == false this simpliy sets the given index, if par == true the index has to be converted to the split index.
	 */
	CUDA_HOST_DEVICE inline void setIndexFromNonParitySplitOrder( lat_index_t latticeIndex )
	{
		if( par == NO_SPLIT )
		{
			index = latticeIndex;
		}
		else if( par == FULL_SPLIT || par == TIMESLICE_SPLIT )
		{
			lat_coord_t site[4];
			lat_index_t parity = 0;
			for( lat_dim_t i = Nd-1; i >= 0; i-- )
			{
				site[i] = latticeIndex % size[i];
				parity += site[i];
				latticeIndex /= size[i];
			}
			index = getIndex( site );
		}
		else
		{
			assert(false);
		}
	}

	/**
	 * Returns the lattice size.
	 * @return lattice size
	 */
	CUDA_HOST_DEVICE inline lat_index_t getSize() const
	{
		return latticeSize;
	}


	/**
	 * TODO
	 * Returns the index within a timeslice
	 * @return lattice index in timeslice
	 */
	CUDA_HOST_DEVICE inline lat_index_t getIndexTimeslice() const
	{
		if( par == FULL_SPLIT )
		{
			assert(false); // TODO
			return -1;
		}
		else if ( par == TIMESLICE_SPLIT )
		{
			return index % getSizeTimeslice();
		}
		else // NO_SPLIT
		{
			return index % getSizeTimeslice();
		}
	}

	/**
	 * Returns the size of a timeslice, i.e. latticeSize/size[0]
	 * @return size of timeslice
	 */
	CUDA_HOST_DEVICE inline lat_index_t getSizeTimeslice() const
	{
		return latticeSize/size[0];
	}

	/**
	 * Returns the size of the lattice in direction i
	 * @return size of timeslice
	 */
	CUDA_HOST_DEVICE inline lat_coord_t getSizeDirection( lat_dim_t i )
	{
		return size[i];
	}

	/**
	 * Sets the index to the neighbour given by mu = 0..(Nd-1) and direction "up": positive direction (up==true), negative direction (up==false)
	 * using a neighbour table stored in memory. On the fly calculation reduces memory traffic but probably needs more registers.
	 */
	CUDA_HOST_DEVICE inline void setNeighbour( lat_dim_t mu, bool up )
	{
		index = nn[(2*mu+up)*getSize()+index];
	}

	/**
	 * Sets the index to the neighbour given by mu = 0..(Nd-1) and direction "up": positive direction (up==true), negative direction (up==false)
	 * using a neighbour table stored in memory. On the fly calculation reduces memory traffic but probably needs more registers.
	 */
	CUDA_HOST_DEVICE inline SiteIndex<Nd, par> getNeighbour( lat_dim_t mu, bool up )
	{
		SiteIndex<Nd, par> site(*this);
		site.setNeighbour( mu, up );
		return site;
	}

	/**
	 * Calculates the table of neighbours. No effort is made to do this in a fast way since this method is invoked only once...
	 * TODO the parameter "nn" should be placed in the constructor.
	 */
	CUDA_HOST_DEVICE inline void calculateNeighbourTable( lat_index_t* nn )
	{
		for( lat_index_t i = 0; i < getSize(); i++ )
		{
			lat_coord_t site[Nd];
			lat_index_t latticeIndex = i;
			for( lat_dim_t j = Nd-1; j >= 0; j-- ) // calculate site vector
			{
				site[j] = latticeIndex % size[j];
				latticeIndex /= size[j];
			}

			lat_coord_t copySite[Nd];
			for( lat_dim_t j = 0; j < Nd; j++ )
			{
				for( int k = 0; k < Nd; k++ )
					copySite[k] = site[k];
				copySite[j]--;
				if( copySite[j] < 0 ) copySite[j] += size[j];
				nn[(2*j)*getSize()+getIndex(site)] = getIndex(copySite);

				for( int k = 0; k < Nd; k++ )
					copySite[k] = site[k];
				copySite[j]++;
				if( copySite[j] >= size[j] ) copySite[j] -= size[j];
				nn[(2*j+1)*getSize()+getIndex(site)] = getIndex(copySite);
			}
		}
	}

	/**
	 * Helper function for calculating the neighbour table.
	 */
	CUDA_HOST_DEVICE inline lat_index_t getIndex( const lat_coord_t site[Nd] )
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
};

}


#endif /* SITEINDEX_HXX_ */
