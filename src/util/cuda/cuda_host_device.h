/*
 * cuda_host_device.h
 *
 *  Created on: Apr 12, 2012
 *      Author: vogt
 */

#ifndef CUDA_HOST_DEVICE_H_
#define CUDA_HOST_DEVICE_H_


#ifdef __CUDA_ARCH__
#define CUDA_HOST_DEVICE __device__ __host__
#define CUDA
#else
#define CUDA_HOST_DEVICE
#endif


// this is quatsch oder?
#ifdef __CUDA_ARCH__
#define CUDA_DEVICE __device__
#else
#define CUDA_DEVICE
#endif

#endif /* CUDA_HOST_DEVICE_H_ */
