/*
 *
 */
#include "DeviceProperties.h"
namespace culgt
{
	int DeviceProperties::getComputeCapability()
	{
		if( !isAvailable )
		{
			update();
		}
		return computeCapability;
	}
	int DeviceProperties::getMaxGridSize()
	{
		if( !isAvailable )
		{
			update();
		}
		return maxGridSize;
	}
	int DeviceProperties::getMaxBlockSize()
	{
		if( !isAvailable )
		{
			update();
		}
		return maxBlockSize;
	}
	std::string DeviceProperties::getName()
	{
		if( !isAvailable )
		{
			update();
		}
		return name;
	}
	int DeviceProperties::getDeviceNumber()
	{
		if( !isAvailable )
		{
			update();
		}
		return deviceNumber;
	}
	void DeviceProperties::update()
	{
		cudaDeviceProp deviceProp;
		cudaGetDevice( &deviceNumber );
		cudaGetDeviceProperties(&deviceProp, deviceNumber );

		computeCapability = deviceProp.major*100+deviceProp.minor*10;
		maxGridSize = deviceProp.maxGridSize[0];
		maxBlockSize = deviceProp.maxThreadsPerBlock;
		name = deviceProp.name;

		isAvailable = true;
	}
	bool DeviceProperties::isAvailable = false;
	int DeviceProperties::computeCapability;
	int DeviceProperties::maxGridSize;
	int DeviceProperties::maxBlockSize;
	std::string DeviceProperties::name;
	int DeviceProperties::deviceNumber;
}

