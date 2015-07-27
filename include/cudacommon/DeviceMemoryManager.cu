/*
 *
 * This is a very basic class that handles device memory allocation. Should be upgrade to allow
 *  - multiple devices
 *
 *  Created on: Mar 6, 2015
 *      Author: vogt
 */

#include "DeviceMemoryManager.h"


using std::ios;

namespace culgt
{

__global__ void clearKernel( void* ptr, size_t size )
{
	for( size_t i = 0; i < size/sizeof(char); ++i )
	{
		((char*)ptr)[i] = '0';
	}
}

bad_alloc_cuda::bad_alloc_cuda( std::string str ) throw() :str(str) { }
bad_alloc_cuda::~bad_alloc_cuda() throw(){};
const char* bad_alloc_cuda::what() const throw()
{
	return str.c_str();
}

MemoryStruct::MemoryStruct( size_t size, std::string info ): size(size), info(info)
{
}
//MemoryStruct::MemoryStruct( const MemoryStruct& m ): size(m.getSize()), info(m.getInfo())
//{
//}

std::string MemoryStruct::getInfo() const
{
	return info;
}

size_t MemoryStruct::getSize() const
{
	return size;
}

void DeviceMemoryManager::registerMalloc( void* pointer, size_t size, std::string description )
{
	info.insert( std::pair<void*, MemoryStruct>( pointer, MemoryStruct( size, description) ) );
	if( verbose )
	{
		ios init(NULL);
		init.copyfmt(std::cout);
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
		std::cout.copyfmt(init);
	}
}

void DeviceMemoryManager::registerFree( void* pointer )
{
	ios init(NULL);
	init.copyfmt(std::cout);
	if( verbose )
	{
		std::cout << "-----------------------------------------------------\n";
		std::cout << "Freeing memory allocated for " << info.at(pointer).getInfo() << "\n";
		std::cout << "Freeing:          " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
				<< (double)info.at(pointer).getSize()/1024./1024. << " MB\n";
	}
	info.erase( pointer );
	if(verbose)
	{
		std::cout << "Total allocated:  " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getMemoryUsageMB() << " MB\n";
		std::cout << "Available memory: " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFreeMemoryMB() << " MB\n";
	}
	std::cout.copyfmt(init);
}

size_t DeviceMemoryManager::getMemoryUsage()
{
	size_t allocatedMemory = 0;
	for ( std::map<void*, MemoryStruct>::iterator it=info.begin(); it!=info.end(); ++it)
	{
		allocatedMemory += it->second.getSize();
	}
   return allocatedMemory;
}

double DeviceMemoryManager::getMemoryUsageMB()
{
	return (double)getMemoryUsage()/1024./1024.;
}

double DeviceMemoryManager::getUnregisteredMemoryMB()
{
	return (double)(getTotalAllocatedMemory()-getMemoryUsage())/1024/1024;
}

size_t DeviceMemoryManager::getFreeMemory()
{
	size_t freeMem;
	size_t total;
	cudaMemGetInfo( &freeMem, &total );
	return freeMem;
}

double DeviceMemoryManager::getFreeMemoryMB()
{
	return (double)getFreeMemory()/1024./1024.;
}

void DeviceMemoryManager::setVerbose()
{
	verbose = true;
}

void DeviceMemoryManager::unsetVerbose()
{
	verbose = false;
}

size_t DeviceMemoryManager::getTotalAllocatedMemory()
{
	size_t freeMem;
	size_t total;
	cudaMemGetInfo( &freeMem, &total );
	return (total-freeMem);
}

size_t DeviceMemoryManager::allocatedMemory = 0;
bool DeviceMemoryManager::verbose = false;
std::map<void*, MemoryStruct> DeviceMemoryManager::info;

}
