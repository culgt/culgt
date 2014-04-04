/**
 * KernelSetup.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef KERNELSETUP_H_
#define KERNELSETUP_H_
#include "LatticeDimension.h"

namespace culgt
{

#define VERIFY_LATTICE_SIZE( dim, index )\
	if( index >= dim.getSize() ) return

template<int Ndim> class KernelSetup
{
public:
	KernelSetup( LatticeDimension<Ndim> dim, bool paritySplit = false, int blockSize = 32 ) : dim( dim ), blockSize(blockSize)
	{
		lat_index_t latSize = dim.getSize();
		if( paritySplit ) latSize /= 2;
		int mod = latSize%blockSize;
		gridSize = latSize/blockSize;
		if( mod > 0 ) gridSize++;
	}
	int getBlockSize() const
	{
		return blockSize;
	}

	int getGridSize() const
	{
		return gridSize;
	}

//	static void verifyLattice( LatticeDimension<Ndim> dim, int index )
//	{
//		if( dim.getLatticeSize() <= index ) return;
//	}

private:
	LatticeDimension<Ndim> dim;
	int gridSize;
	int blockSize;
};

}
#endif /* KERNELSETUP_H_ */
