/**
 * KernelSetup.h
 *
 *
 *	workSize is the number of threads that are needed per working block. KernelSetup tries to use the smallest multiple of workSize threads
 *	that are compatible with DeviceProperties::getMaxGridSize().
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef KERNELSETUP_H_
#define KERNELSETUP_H_
#include "LatticeDimension.h"
#include <iostream>
#include <string>
#include <sstream>
#include "../cudacommon/DeviceProperties.h"

using std::stringstream;

namespace culgt
{

#define VERIFY_LATTICE_SIZE( dim, index )\
	if( index >= dim.getSize() ) return

template<int Ndim> class KernelSetup
{
public:
	KernelSetup( LatticeDimension<Ndim> dim, bool paritySplit = false, int minSitesPerBlock = 32, int threadsPerSite = 1 ) : dim( dim )
	{
		int smallestSitesPerBlock = getSmallestBlockSize( paritySplit, minSitesPerBlock, threadsPerSite );
		if( smallestSitesPerBlock > 0 )
			initSetup( paritySplit, smallestSitesPerBlock, threadsPerSite );
		else
		{
			// TODO ugly style: print message and let the kernel call die...
			std::cout << "No suitable KernelSetup found (problem too large). Need to upgrade KernelSetup to deal with that..." << std::endl;
			this->gridSize = -1;
			this->blockSize = -1;
		}

	}
	int getBlockSize() const
	{
		return blockSize;
	}

	int getGridSize() const
	{
		return gridSize;
	}

private:
	LatticeDimension<Ndim> dim;
	int gridSize;
	int blockSize;

	bool initSetup( bool paritySplit, int smallestSitesPerBlock, int threadsPerSite )
	{
		this->blockSize = smallestSitesPerBlock*threadsPerSite;
		lat_index_t latSize = dim.getSize();
		if( paritySplit ) latSize /= 2;
		int mod = latSize%smallestSitesPerBlock;
		gridSize = latSize/smallestSitesPerBlock;
		if( mod > 0 ) gridSize++;
		if( blockSize <= DeviceProperties::getMaxBlockSize() && gridSize <= DeviceProperties::getMaxGridSize() ) return true;
		else return false;
	}

	int getSmallestBlockSize( bool paritySplit, int minSitesPerBlock, int threadsPerBlock )
	{
		lat_index_t latSize = dim.getSize();
		if( paritySplit ) latSize /= 2;

		int sitesPerBlock = minSitesPerBlock;
		while( true )
		{
			if( latSize/sitesPerBlock < DeviceProperties::getMaxGridSize() - 1 )
			{
				break;
			}
			else
			{
				sitesPerBlock += minSitesPerBlock;
			}
		}
		if( sitesPerBlock*threadsPerBlock <= DeviceProperties::getMaxBlockSize() ) return sitesPerBlock;
		else return -1;
	}

};

}
#endif /* KERNELSETUP_H_ */
