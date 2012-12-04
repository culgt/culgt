/**
 * Defines global constants in host and device memory.
 * They are defined in different namespaces for host and device respectively.
 *
 * TODO: use a class to hide the pointers and enforce use of the getters.
 */

#ifndef GLOBALCONSTANTS_HXX_
#define GLOBALCONSTANTS_HXX_

#include "../lattice/cuda/cuda_host_device.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include <cuda.h>

#ifdef _X_
const lat_coord_t Nx = _X_;
#else
#error "Define X (the lattice size in x-direction)"
#endif
#ifdef _Y_
const lat_coord_t Ny = _Y_;
#else
const lat_coord_t Ny = _X_;
#warning "Size in y-direction not set. Using size of x."
#endif
#ifdef _Z_
const lat_coord_t Nz = _Z_;
#else
const lat_coord_t Nz = _X_;
#warning "Size in z-direction not set. Using size of x."
#endif
#ifdef _T_
const lat_coord_t Nt = _T_;
#else
#error "Define T (the lattice size in t-direction)"
#endif

// number of sites per block
#ifdef _NSB_
const int NSB = _NSB_;
const int NSB2= 2*NSB;
const int NSB3= 3*NSB;
const int NSB4= 4*NSB;
#else
const int NSB = 32;
const int NSB2= 2*NSB;
const int NSB3= 3*NSB;
const int NSB4= 4*NSB;
//#error "Define NSB (number of lattice sites per ThreadBlock (e.g. 32))"
#endif

/**
 * When using SIZE/SIZE_TIMESLICE in kernels: look at your register spilling. In a test the spilling is much increased when using this instead of locally defining the size array!
 */
namespace DEVICE_CONSTANTS
{
//	extern __constant__ lat_coord_t SIZE[4];
//	extern __constant__ lat_coord_t SIZE_TIMESLICE[4];
	__constant__ lat_coord_t SIZE[4]  = {Nt,Nx,Ny,Nz};
	__constant__ lat_coord_t SIZE_TIMESLICE[4] = {1,Nx,Ny,Nz};
}


namespace HOST_CONSTANTS
{
	const lat_coord_t SIZE[4]  = {Nt,Nx,Ny,Nz};
	const lat_coord_t SIZE_TIMESLICE[4]  = {1,Nx,Ny,Nz};

	extern lat_coord_t* PTR_TO_DEVICE_SIZE;
	extern lat_coord_t* PTR_TO_DEVICE_SIZE_TIMESLICE;
	extern bool isInitialized;

	void init();

	lat_coord_t* getPtrToDeviceSize();
	lat_coord_t* getPtrToDeviceSizeTimeslice();
}


// TODO the following is what I want to do, but seems like not supported on CUDA (or bug? or error on my side?)

//struct DEVICE_CONSTANTS
//{
//	static __constant__ lat_coord_t SIZE[];
//};
//
//lat_coord_t DEVICE_CONSTANTS::SIZE[]  = {Nt,Nx,Ny,Nz};

//struct HOST_CONSTANTS
//{
//	static const lat_coord_t SIZE[];
//};
//
//const lat_coord_t HOST_CONSTANTS::SIZE[] = {Nt,Nx,Ny,Nz};


#endif /* GLOBALCONSTANTS_HXX_ */
