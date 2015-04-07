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

namespace culgt
{

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
//
//size_t DeviceMemoryManager::allocatedMemory = 0;
//bool DeviceMemoryManager::verbose = true;
//std::map<void*, MemoryStruct> DeviceMemoryManager::info;

}

#endif /* DEVICEMEMORYMANAGER_H_ */