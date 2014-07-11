/**
 * DeviceReader.h
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */

#ifndef DEVICECOMMUNICATOR_H_
#define DEVICECOMMUNICATOR_H_

#include "cuda_error.h"

namespace culgt
{

template<typename T> class DeviceCommunicator
{
public:
	static T getValue( T* devPtr, int index )
	{
		T hostVal;// = (T)0;
		CUDA_SAFE_CALL( cudaMemcpy( &hostVal, &devPtr[index], sizeof(T), cudaMemcpyDeviceToHost ), "DeviceCommunicator memCpy to host" );
		return hostVal;
	}

	static void setValue( T* devPtr, int index, T val )
	{
		CUDA_SAFE_CALL( cudaMemcpy(  &devPtr[index], &val, sizeof(T), cudaMemcpyHostToDevice ), "DeviceCommunicator memCpy to device" );
	}
};

}

#endif /* DEVICECOMMUNICATOR_H_ */
