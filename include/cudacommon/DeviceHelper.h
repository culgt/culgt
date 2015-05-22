/*
 * DeviceHelper.h
 *
 *  Created on: May 22, 2015
 *      Author: vogt
 */

#ifndef DEVICEHELPER_H_
#define DEVICEHELPER_H_

class DeviceHelper
{
public:
	static int selectAvailableDevice()
	{
		int nDevices;
		cudaGetDeviceCount( &nDevices );

		for( int i = 0; i < nDevices; ++i )
		{
			cudaError_t error = cudaSetDevice( i );
			if( error != cudaSuccess ) // I would expect that this should already be != success if device is busy in exclusive mode but it did not work
				continue;
			cudaFree(0);
			cudaDeviceSynchronize();
			error = cudaGetLastError();
			if( error == cudaSuccess )
					return i;
		}
		return -1;
	}
};


#endif /* DEVICEHELPER_H_ */
