/**
 * Lattice Site using a 1 dimensional index.
 *
 * This implementation of "Site" uses the 1 dimensional index to store the current site. It is favorable in CUDA implementations
 * since it uses fewer variables (exactly one) to store the site compared to "SiteCoord" which uses Ndim variables.
 * This results in fewer registers which is very important for CUDA optimization.
 * Some functions are more costly in this implementation, like operator[] (access to coordinates), other functions are much faster,
 * like getLatticeIndex() (which simply return the member variable).
 *
 * Possible Optimizations:
 *  - I don't know if the compiler is smart enough to place size[] in constant memory. Check this and think about a implementation
 *    where size is stored in constant memory. (For example pass a constant array by reference.)
 *  - Think about the Quandt style neighbour table.
 *  - Think about calculating neighbour on the fly.
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-05-02
 */

#ifndef SITEINDEX_HXX_
#define SITEINDEX_HXX_

#include "cuda/cuda_host_device.h"
#include "datatype/lattice_typedefs.h"
#include <assert.h>

/**
 * The template parameter "bool par" defines normal indexing (par==false) or parity split indexing (par==true).
 * Parity split indexing means that the first half of the indices are the even sites ((x+y+z+...)%2 == 0) and the second half the odd sites.
 * It has to be used for CUDA implementations to ensure coalesced memory reads for neighbouring (concerning its index) array elements.
 *
 * @param template Nd: Dimension of the lattice
 * @param template par
 */
template<lat_dim_t Nd, ParityType par> class SiteIndex
{
public:
	CUDA_HOST_DEVICE inline SiteIndex( const lat_coord_t size[Nd] );
	CUDA_HOST_DEVICE inline SiteIndex( const SiteIndex<Nd,par> &s);
	CUDA_HOST_DEVICE inline virtual ~SiteIndex();
	CUDA_HOST_DEVICE inline lat_coord_t operator[](lat_dim_t i);
	CUDA_HOST_DEVICE inline lat_index_t getLatticeIndex();
	CUDA_HOST_DEVICE inline void setLatticeIndex( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline void setLatticeIndexTimeslice( lat_index_t latticeIndex, lat_coord_t t );
	CUDA_HOST_DEVICE inline void setLatticeIndexFromParitySplitOrder( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline void setLatticeIndexFromNonParitySplitOrder( lat_index_t latticeIndex );
	CUDA_HOST_DEVICE inline lat_index_t getLatticeSize();
	CUDA_HOST_DEVICE inline lat_index_t getLatticeIndexTimeslice();
	CUDA_HOST_DEVICE inline lat_index_t getLatticeSizeTimeslice();
	CUDA_HOST_DEVICE inline lat_coord_t getLatticeSizeDirection( lat_dim_t i );
	CUDA_HOST_DEVICE inline void setNeighbour( lat_dim_t direction, bool up );

	CUDA_HOST_DEVICE inline void calculateNeighbourTable( lat_index_t* nn );

	static const lat_dim_t Ndim = Nd;
	lat_index_t* nn;

	lat_coord_t size[Nd];
	CUDA_HOST_DEVICE inline lat_index_t getLatticeIndex( const lat_coord_t site[Nd] );
private:
	lat_index_t index;
	lat_index_t latticeSize;
};


template <lat_dim_t Nd, ParityType par> SiteIndex<Nd, par>::SiteIndex( const lat_coord_t size[Nd] )
{
	latticeSize = 1;
	for( lat_dim_t i = 0; i < Nd; i++ )
	{
		this->size[i] = size[i];
		latticeSize *= size[i];
	}
}

template <lat_dim_t Nd, ParityType par> SiteIndex<Nd, par>::SiteIndex( const SiteIndex<Nd,par> &s )
{
	this->index = s.index;
	latticeSize = 1;
	for( lat_dim_t i = 0; i < Nd; i++ )
	{
		this->size[i] = s.size[i];
		latticeSize *= size[i];
	}
	this->nn = s.nn;
}

template <lat_dim_t Nd, ParityType par> SiteIndex<Nd, par>::~SiteIndex()
{
}

/**
 * TODO
 * Returns the i-coordinate of the current site.
 */
template<lat_dim_t Nd, ParityType par> lat_coord_t SiteIndex<Nd,par>::operator[](lat_dim_t i)
{
	lat_index_t temp = index;

	if( par == FULL_SPLIT )
	{
		assert(false); // TODO
		return -1;
	}
	else if ( par == TIMESLICE_SPLIT )
	{
		assert(false); // TODO
		return -1;
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
template<lat_dim_t Nd, ParityType par> lat_index_t SiteIndex<Nd, par>::getLatticeIndex()
{
	return index;
}

/**
 * Sets the lattice index.
 * @param lattice index
 * @return void
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::setLatticeIndex( lat_index_t latticeIndex )
{
	index = latticeIndex;
}

/**
 * Sets the lattice index.
 * @param lattice index
 * @return void
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::setLatticeIndexTimeslice( lat_index_t latticeIndex, lat_coord_t t )
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
		index = t*getLatticeSizeTimeslice() + latticeIndex;
	}
}

/**
 * // TODO
 * If par == true this simpliy sets the given index, if par == false the index has to be converted to the non-split index.
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::setLatticeIndexFromParitySplitOrder( lat_index_t latticeIndex )
{
	// TODO
	assert(false);
}

/**
 * // TODO
 * If par == false this simpliy sets the given index, if par == true the index has to be converted to the split index.
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::setLatticeIndexFromNonParitySplitOrder( lat_index_t latticeIndex )
{
	// TODO
	assert(false);
}

/**
 * Returns the lattice size.
 * @return lattice size
 */
template<lat_dim_t Nd, ParityType par> lat_index_t SiteIndex<Nd, par>::getLatticeSize()
{
	return latticeSize;
}


/**
 * TODO
 * Returns the index within a timeslice
 * @return lattice index in timeslice
 */
template<lat_dim_t Nd, ParityType par> lat_index_t SiteIndex<Nd, par>::getLatticeIndexTimeslice()
{
	if( par == FULL_SPLIT )
	{
		assert(false); // TODO
		return -1;
	}
	else if ( par == TIMESLICE_SPLIT )
	{
		assert(false); // TODO
		return -1;
	}
	else // NO_SPLIT
	{
		return index % getLatticeSizeTimeslice();
	}
}

/**
 * Returns the size of a timeslice, i.e. latticeSize/size[0]
 * @return size of timeslice
 */
template<lat_dim_t Nd, ParityType par> lat_index_t SiteIndex<Nd, par>::getLatticeSizeTimeslice()
{
	return latticeSize/size[0];
}

/**
 * Returns the size of the lattice in direction i
 * @return size of timeslice
 */
template<lat_dim_t Nd, ParityType par> lat_coord_t SiteIndex<Nd, par>::getLatticeSizeDirection( lat_dim_t i )
{
	return size[i];
}

/**
 * Sets the index to the neighbour given by mu = 0..(Nd-1) and direction "up": positive direction (up==true), negative direction (up==false)
 * using a neighbour table stored in memory. On the fly calculation reduces memory traffic but probably needs more registers.
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::setNeighbour( lat_dim_t mu, bool up )
{
	index = nn[(2*mu+up)*getLatticeSize()+index];
}

/**
 * Calculates the table of neighbours. No effort is made to do this in a fast way since this method is invoked only once...
 * TODO the parameter "nn" should be placed in the constructor.
 */
template<lat_dim_t Nd, ParityType par> void SiteIndex<Nd, par>::calculateNeighbourTable( lat_index_t* nn )
{
	for( lat_index_t i = 0; i < getLatticeSize(); i++ )
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
			nn[(2*j)*getLatticeSize()+getLatticeIndex(site)] = getLatticeIndex(copySite);

			for( int k = 0; k < Nd; k++ )
				copySite[k] = site[k];
			copySite[j]++;
			if( copySite[j] >= size[j] ) copySite[j] -= size[j];
			nn[(2*j+1)*getLatticeSize()+getLatticeIndex(site)] = getLatticeIndex(copySite);
		}
	}
}

/**
 * Helper function for calculating the neighbour table.
 */
template<lat_dim_t Nd, ParityType par> lat_index_t SiteIndex<Nd, par>::getLatticeIndex( const lat_coord_t site[Nd] )
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
			return index / 2 + getLatticeSize()/2;
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
			index = index / 2 + getLatticeSizeTimeslice()/2;
		}
		return index + site[0] * getLatticeSizeTimeslice();

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




#endif /* SITEINDEX_HXX_ */
