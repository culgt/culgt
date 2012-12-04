/**
 * Defines global constants in host and device memory.
 * They are defined in different namespaces for host and device respectively.
 *
 * TODO: use a class to hide the pointers and enforce use of the getters.
 */


#include "GlobalConstants.h"

lat_coord_t* HOST_CONSTANTS::PTR_TO_DEVICE_SIZE;
lat_coord_t* HOST_CONSTANTS::PTR_TO_DEVICE_SIZE_TIMESLICE;
bool HOST_CONSTANTS::isInitialized = false;
//
//__constant__ lat_coord_t DEVICE_CONSTANTS::SIZE[4]  = {Nt,Nx,Ny,Nz};
//__constant__ lat_coord_t DEVICE_CONSTANTS::SIZE_TIMESLICE[4] = {1,Nx,Ny,Nz};

void HOST_CONSTANTS::init()
{
	cudaGetSymbolAddress((void **)&PTR_TO_DEVICE_SIZE, DEVICE_CONSTANTS::SIZE );
	cudaGetSymbolAddress((void **)&PTR_TO_DEVICE_SIZE_TIMESLICE, DEVICE_CONSTANTS::SIZE_TIMESLICE );
	isInitialized = true;
}

lat_coord_t* HOST_CONSTANTS::getPtrToDeviceSize()
{
	if( isInitialized ) return PTR_TO_DEVICE_SIZE;
	else
	{
		init();
		return PTR_TO_DEVICE_SIZE;
	}
}

lat_coord_t* HOST_CONSTANTS::getPtrToDeviceSizeTimeslice()
{
	if( isInitialized ) return PTR_TO_DEVICE_SIZE_TIMESLICE;
	else
	{
		init();
		return PTR_TO_DEVICE_SIZE_TIMESLICE;
	}
}

