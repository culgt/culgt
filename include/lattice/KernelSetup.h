/**
 * KernelSetup.h
 *
 *
 *	workSize is the number of threads that are needed per working block. KernelSetup tries to use the smallest multiple of workSize threads
 *	that are compatible with MAX_GRIDSIZE.
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
// TODO This are the FERMI settings... Should deal somewhere with deviceProperties...
int MAX_GRIDSIZE = DeviceProperties::getMaxGridSize();
int MAX_BLOCKSIZE = DeviceProperties::getMaxBlockSize();

//class KernelSetupException: public std::exception
//{
//public:
//	~KernelSetupException() throw() {};
//	KernelSetupException( std::string msg ) : msg(msg){};
//	virtual const char* what() const throw()
//	{
//		return msg.c_str();
//	}
//private:
//	std::string msg;
//};


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


//		if( blockSize == -1 )
//		{
//			int smallestBlockSize = getSmallestBlockSize( paritySplit );
//			if( smallestBlockSize > 0 )
//				initSetup( paritySplit, smallestBlockSize );
//			else
//			{
//				// TODO ugly style: print message and let the kernel call die...
//				std::cout << "No suitable KernelSetup found (problem too large). Need to upgrade KernelSetup to deal with that..." << std::endl;
//				this->gridSize = -1;
//				this->blockSize = -1;
//
//
////				throw KernelSetupException( "No suitable KernelSetup found (problem too large). Need to upgrade KernelSetup to deal with that..." );
//			}
//		}
//		else
//		{
//			if( !initSetup( paritySplit, blockSize ) )
//			{
//				// TODO ugly style: print message and let the kernel call die...
//				std::cout << "No suitable KernelSetup possible with BlockSize = "  << blockSize << std::endl;
//				this->gridSize = -1;
//				this->blockSize = -1;
//
////				stringstream temp;
////				temp << "No suitable KernelSetup possible with BlockSize = ";
////				temp << blockSize;
////				throw KernelSetupException( temp.str() );
//			}
//		}

//		std::cout << "KernelSetup: " << this->gridSize << "/" << this->blockSize << std::endl;
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
		if( blockSize <= MAX_BLOCKSIZE && gridSize <= MAX_GRIDSIZE ) return true;
		else return false;
	}

	int getSmallestBlockSize( bool paritySplit, int minSitesPerBlock, int threadsPerBlock )
	{
		lat_index_t latSize = dim.getSize();
		if( paritySplit ) latSize /= 2;

		int sitesPerBlock = minSitesPerBlock;
		while( true )
		{
			if( latSize/sitesPerBlock < MAX_GRIDSIZE - 1 )
			{
				break;
			}
			else
			{
				sitesPerBlock += minSitesPerBlock;
			}
		}
		if( sitesPerBlock*threadsPerBlock <= MAX_BLOCKSIZE ) return sitesPerBlock;
		else return -1;
	}

};

}
#endif /* KERNELSETUP_H_ */
