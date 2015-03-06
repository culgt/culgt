/*
 * DeviceMemoryManager.h
 *
 * This is a very basic class that handles device memory allocation. Should be upgrade to allow
 *  - multiple devices
 *  - storing information in a map for each allocation
 *    (this allows for freeing with correct size by just giving the pointer and keep a description of what the memory is used for...)
 *
 *  Created on: Mar 6, 2015
 *      Author: vogt
 */

#ifndef DEVICEMEMORYMANAGER_H_
#define DEVICEMEMORYMANAGER_H_
#include <iomanip>

class DeviceMemoryManager
{
public:
	template<typename T> static void malloc( T** pointerToPointer, size_t size )
	{
		CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, size ) , "cudaMalloc in DeviceMemoryManager" );
		registerMalloc( size );
	}

	template<typename T> static void free( T* pointerToPointer, size_t size )
	{
		CUDA_SAFE_CALL( cudaFree( pointerToPointer ) , "cudaFree in DeviceMemoryManager" );
		registerFree( size );
	}

	static void registerMalloc( size_t allocatedMemory )
	{
		DeviceMemoryManager::allocatedMemory += allocatedMemory;
		if( verbose )
		{
			std::cout << "Allocating:       " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
					<< (double)allocatedMemory/1024./1024. << " MB\n";
			std::cout << "Total allocated:  " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< (double)DeviceMemoryManager::allocatedMemory/1024./1024. << " MB\n";
			std::cout << "Available memory: " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getFreeMemoryMB() << " MB\n";
		}
	}
	static void registerFree( size_t freedMemory )
	{
		DeviceMemoryManager::allocatedMemory -= freedMemory;
		if( verbose )
		{
			std::cout << "Freeing:          " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
					<< (double)freedMemory/1024./1024. << " MB\n";
			std::cout << "Total allocated:  " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< (double)DeviceMemoryManager::allocatedMemory/1024./1024. << " MB\n";
			std::cout << "Available memory: " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getFreeMemoryMB() << " MB\n";
		}
	}
	static double getMemoryUsageMB()
	{
		return (double)allocatedMemory/1024./1024.;
	}
	static double getFreeMemoryMB()
	{
		size_t freeMom;
		size_t total;
		cudaMemGetInfo( &freeMom, &total );
		return (double)freeMom/1024./1024.;
	}
private:
	static size_t allocatedMemory;
	static bool verbose;
};

size_t DeviceMemoryManager::allocatedMemory = 0;
bool DeviceMemoryManager::verbose = true;

#endif /* DEVICEMEMORYMANAGER_H_ */
