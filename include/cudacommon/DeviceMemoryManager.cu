/*
 *
 * This is a very basic class that handles device memory allocation. Should be upgrade to allow
 *  - multiple devices
 *
 *  Created on: Mar 6, 2015
 *      Author: vogt
 */

#include "DeviceMemoryManager.h"

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
		std::ios init(NULL);
		init.copyfmt(std::cout);
		std::cout << "-----------------------------------------------------\n";
		std::cout << "Allocating memory for " << description << "\n";
		std::cout << "Allocating:         " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( size ) << "\n";
		std::cout << "Total allocated:    " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( getMemoryUsage() ) << "\n";
		std::cout << "Unregistered memory " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( getUnregisteredMemory() ) << "\n";
		std::cout << "Available memory:     " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( getFreeMemory() ) << "\n";
		std::cout.copyfmt(init);
	}
}

std::string DeviceMemoryManager::makeFormattedBytesString(double size, std::string unit)
{
	std::ostringstream out;
	out << size << " " << unit;
	return out.str();
}

std::string DeviceMemoryManager::getFormattedBytes( size_t size )
{
	if( size > 10*GIGABYTE )
	{
		return makeFormattedBytesString( (double)size/(double)GIGABYTE, "GB" );
	}
	else if( size > 10*MEGABYTE )
	{
		return makeFormattedBytesString( (double)size/(double)MEGABYTE, "MB" );
	}
	else if( size > 10*KILOBYTE )
	{
		return makeFormattedBytesString( (double)size/(double)KILOBYTE, "KB" );
	}
	else
	{
		return makeFormattedBytesString( (double)size, "Bytes" );
	}

}

void DeviceMemoryManager::registerFree( void* pointer )
{
	std::ios init(NULL);
	init.copyfmt(std::cout);
	if( verbose )
	{
		std::cout << "-----------------------------------------------------\n";
		std::cout << "Freeing memory allocated for " << info.at(pointer).getInfo() << "\n";
		std::cout << "Freeing:          " << std::fixed << std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( info.at(pointer).getSize() ) << "\n";
	}
	info.erase( pointer );
	if(verbose)
	{
		std::cout << "Total allocated:  " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( getMemoryUsage() ) << "\n";
		std::cout << "Available memory: " <<  std::setw( 15 ) <<  std::setprecision(3)
				<< getFormattedBytes( getFreeMemory() ) << " \n";
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

size_t DeviceMemoryManager::getUnregisteredMemory()
{
	return getTotalAllocatedMemory()-getMemoryUsage();
}

size_t DeviceMemoryManager::getFreeMemory()
{
	size_t freeMem;
	size_t total;
	cudaMemGetInfo( &freeMem, &total );
	return freeMem;
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

