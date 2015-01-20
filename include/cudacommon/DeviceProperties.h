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
	static int getComputeCapability();
	static int getMaxGridSize();
	static int getMaxBlockSize();
	static std::string getName();
	static int getDeviceNumber();
private:
	static void update();

	static bool isAvailable;
	static int computeCapability;
	static int maxGridSize;
	static int maxBlockSize;
	static std::string name;
	static int deviceNumber;
};

}


#endif /* DEVICEPROPERTIES_H_ */
