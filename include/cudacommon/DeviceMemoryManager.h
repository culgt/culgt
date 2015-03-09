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
#include <string>
#include <map>

namespace culgt
{

class bad_alloc_cuda : public std::bad_alloc
{
public:
	bad_alloc_cuda( std::string str ) throw() :str(str) { }

	~bad_alloc_cuda() throw(){};

	const char* what() const throw()
	{
		return str.c_str();
	}

private:
	std::string str;
};

class MemoryStruct
{
public:
	MemoryStruct( size_t size, std::string info ): size(size), info(info)
	{
	}

	std::string getInfo() const
	{
		return info;
	}

	size_t getSize() const
	{
		return size;
	}

private:
	size_t size;
	std::string info;
};


class DeviceMemoryManager
{
public:
	template<typename T> static void malloc( T** pointerToPointer, size_t size, std::string description )
	{
		if( size > getFreeMemory() )
		{
			throw bad_alloc_cuda( "out of memory" );
		}
		CUDA_SAFE_CALL( cudaMalloc( pointerToPointer, size ) , "cudaMalloc in DeviceMemoryManager" );
		registerMalloc( (void*)*pointerToPointer, size, description );
	}

	template<typename T> static void free( T* pointer )
	{
		CUDA_SAFE_CALL( cudaFree( pointer ) , "cudaFree in DeviceMemoryManager" );
		registerFree( pointer );
	}

	static void registerMalloc( void* pointer, size_t size, std::string description )
	{
		info.insert( std::pair<void*, MemoryStruct>( pointer, MemoryStruct( size, description) ) );
		if( verbose )
		{
			std::cout << "-----------------------------------------------------\n";
			std::cout << "Allocating memory for " << description << "\n";
			std::cout << "Allocating:         " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
					<< (double)size/1024./1024. << " MB\n";
			std::cout << "Total allocated:    " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getMemoryUsageMB() << " MB\n";
			std::cout << "Unregistered memory " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getUnregisteredMemoryMB() << "MB\n";
			std::cout << "Available memory:     " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getFreeMemoryMB() << " MB\n";
		}
	}

	static void registerFree( void* pointer )
	{
		if( verbose )
		{
			std::cout << "-----------------------------------------------------\n";
			std::cout << "Freeing memory allocated for " << info.at(pointer).getInfo() << "\n";
			std::cout << "Freeing:          " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
					<< (double)info.at(pointer).getSize()/1024./1024. << " MB\n";

			info.erase( pointer );
			std::cout << "Total allocated:  " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getMemoryUsageMB() << " MB\n";
			std::cout << "Available memory: " <<  std::setw( 15 ) <<  std::setprecision(3)
					<< getFreeMemoryMB() << " MB\n";
		}
	}

	static size_t getMemoryUsage()
	{
		size_t allocatedMemory = 0;
		for ( std::map<void*, MemoryStruct>::iterator it=info.begin(); it!=info.end(); ++it)
		{
			allocatedMemory += it->second.getSize();
		}
	   return allocatedMemory;
	}

	static double getMemoryUsageMB()
	{
		return (double)getMemoryUsage()/1024./1024.;
	}

	static double getUnregisteredMemoryMB()
	{
		return (double)(getTotalAllocatedMemory()-getMemoryUsage())/1024/1024;
	}

private:
	static size_t allocatedMemory;
	static bool verbose;
	static std::map<void*, MemoryStruct> info;

	static size_t getFreeMemory()
	{
		size_t freeMem;
		size_t total;
		cudaMemGetInfo( &freeMem, &total );
		return freeMem;
	}

	static double getFreeMemoryMB()
	{
		return (double)getFreeMemory()/1024./1024.;
	}

	static size_t getTotalAllocatedMemory()
	{
		size_t freeMem;
		size_t total;
		cudaMemGetInfo( &freeMem, &total );
		return (total-freeMem);
	}
};

size_t DeviceMemoryManager::allocatedMemory = 0;
bool DeviceMemoryManager::verbose = true;
std::map<void*, MemoryStruct> DeviceMemoryManager::info;

}

#endif /* DEVICEMEMORYMANAGER_H_ */
