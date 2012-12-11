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
