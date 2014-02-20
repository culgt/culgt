/**
 */

#ifndef CUDA_HOST_DEVICE_H_
#define CUDA_HOST_DEVICE_H_


#ifdef __CUDA_ARCH__
#define CUDA_HOST_DEVICE __device__ __host__
#define CUDA
#else
#define CUDA_HOST_DEVICE
#endif

#endif /* CUDA_HOST_DEVICE_H_ */
