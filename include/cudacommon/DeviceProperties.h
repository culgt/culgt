/*
 * DeviceProperties.h
 *
 *  Created on: Jan 14, 2015
 *      Author: vogt
 */

#ifndef DEVICEPROPERTIES_H_
#define DEVICEPROPERTIES_H_

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
private:
	static void update()
	{
		cudaDeviceProp deviceProp;
		int selectedDeviceNumber;
		cudaGetDevice( &selectedDeviceNumber );
		cudaGetDeviceProperties(&deviceProp, selectedDeviceNumber );

		computeCapability = deviceProp.major*100+deviceProp.minor+10;
		maxGridSize = deviceProp.maxGridSize[0];
		maxBlockSize = deviceProp.maxThreadsPerBlock;

		isAvailable = true;
	}
	static bool isAvailable;
	static int computeCapability;
	static int maxGridSize;
	static int maxBlockSize;
};
bool DeviceProperties::isAvailable = false;
int DeviceProperties::computeCapability;
int DeviceProperties::maxGridSize;
int DeviceProperties::maxBlockSize;

}


#endif /* DEVICEPROPERTIES_H_ */
