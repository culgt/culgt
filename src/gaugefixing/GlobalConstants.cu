/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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

