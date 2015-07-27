/*
 * DeviceMemoryManager.h
 *
 * This is a very basic class that handles device memory allocation. Should be upgrade to allow
 *  - multiple devices
 *
 *  Created on: Mar 6, 2015
 *      Author: vogt
 */

#ifndef DEVICEMEMORYMANAGER_H_
#define DEVICEMEMORYMANAGER_H_
#include <iomanip>
#include <string>
#include <map>
#include <iostream>
#include <iosfwd>

#include "cudacommon/cuda_error.h"

namespace culgt
{

__global__ void clearKernel( void* ptr, size_t size );

class bad_alloc_cuda : public std::bad_alloc
{
public:
	bad_alloc_cuda( std::string str ) throw();
	~bad_alloc_cuda() throw();
	const char* what() const throw();
private:
	std::string str;
};

class MemoryStruct
{
public:
	MemoryStruct( size_t size, std::string info );
	std::string getInfo() const;
	size_t getSize() const;
private:
	size_t size;
	std::string info;
};

class DeviceMemoryManager
{
public:
	// need inline to avoid "multiple definitions" when using in seperated compilation
	template<typename T> static inline void malloc( T** pointerToPointer, size_t size, std::string description )
	{
		if( size > getFreeMemory() )
		{
			throw bad_alloc_cuda( "out of memory" );
		}
		CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, size ) , "cudaMalloc in DeviceMemoryManager" );
		registerMalloc( (void*)*pointerToPointer, size, description );
	}

	template<typename T> static inline void free( T* pointer )
	{
		CUDA_SAFE_CALL( cudaFree( pointer ) , "cudaFree in DeviceMemoryManager" );
		registerFree( pointer );
	}

	template<typename T> static inline void clear( T* pointer )
	{
		std::map<void*, MemoryStruct>::iterator it = info.find( pointer );
		clearKernel<<<1,1>>>((void*)pointer, it->second.getSize() );
	}

	static void registerMalloc( void* pointer, size_t size, std::string description );
	static void registerFree( void* pointer );
	static size_t getMemoryUsage();
	static double getMemoryUsageMB();
	static double getUnregisteredMemoryMB();
	static void setVerbose();
	static void unsetVerbose();


	static bool verbose;
private:
	static size_t allocatedMemory;
	static std::map<void*, MemoryStruct> info;
	static size_t getFreeMemory();
	static double getFreeMemoryMB();
	static size_t getTotalAllocatedMemory();
};

}

#endif /* DEVICEMEMORYMANAGER_H_ */
