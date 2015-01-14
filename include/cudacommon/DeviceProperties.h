/*
 * DeviceProperties.h
 *
 *  Created on: Jan 14, 2015
 *      Author: vogt
 */

#ifndef DEVICEPROPERTIES_H_
#define DEVICEPROPERTIES_H_
#include <string>

namespace culgt
{

class DeviceProperties
{
public:
	static int getComputeCapability()
	{
		if( !isAvailable )
		{
			update();
		}
		return computeCapability;
	}
	static int getMaxGridSize()
	{
		if( !isAvailable )
		{
			update();
		}
		return maxGridSize;
	}
	static int getMaxBlockSize()
	{
		if( !isAvailable )
		{
			update();
		}
		return maxBlockSize;
	}
	static std::string getName()
	{
		if( !isAvailable )
		{
			update();
		}
		return name;
	}
	static int getDeviceNumber()
	{
		if( !isAvailable )
		{
			update();
		}
		return deviceNumber;
	}
private:
	static void update()
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
	static bool isAvailable;
	static int computeCapability;
	static int maxGridSize;
	static int maxBlockSize;
	static std::string name;
	static int deviceNumber;
};
bool DeviceProperties::isAvailable = false;
int DeviceProperties::computeCapability;
int DeviceProperties::maxGridSize;
int DeviceProperties::maxBlockSize;
std::string DeviceProperties::name;
int DeviceProperties::deviceNumber;

}


#endif /* DEVICEPROPERTIES_H_ */
