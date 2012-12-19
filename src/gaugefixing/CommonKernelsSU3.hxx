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

#ifndef COMMONKERNELSSU3_HXX_
#define COMMONKERNELSSU3_HXX_

#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "GlobalConstants.h"


// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.
namespace COMKSU3
{
static const int Ndim = 4; // TODO why here?
static const int Nc = 3;
__global__ void projectSU3( Real *U, lat_coord_t* ptrToDeviceSize );
__global__ void setHot( Real *U, lat_coord_t* ptrToDeviceSize, int rngSeed, int rngCounter );
}

class CommonKernelsSU3
{
public:
//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );

	// TODO remove static and make the init in constructor
	static void initCacheConfig()
	{
		cudaFuncSetCacheConfig( COMKSU3::projectSU3, cudaFuncCachePreferL1 );
	}

	static void projectSU3( int a, int b, Real *U, lat_coord_t* ptrToDeviceSize )
	{
		COMKSU3::projectSU3<<<a,b>>>( U, ptrToDeviceSize );
	};
	
	static void setHot( int a, int b, Real *U, lat_coord_t* ptrToDeviceSize, int rngSeed, int rngCounter )
	{
		COMKSU3::setHot<<<a,b>>>( U, ptrToDeviceSize, rngSeed, rngCounter );
	};

private:
};





namespace COMKSU3
{

__global__ void projectSU3( Real *U, lat_coord_t* ptrToDeviceSize )
{
	typedef GpuPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;
	typedef Link<GpuTimeslice,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(ptrToDeviceSize );
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	// TODO maybe speed up by making link local: Test this.
	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		globUp.projectSU3(); // IMPORTANT: Currently this kernel is used for reconstructing third line in the end. Be aware of this when changing something.
	}
}

__global__ void setHot( Real *U, lat_coord_t* ptrToDeviceSize, int rngSeed, int rngCounter )
{
	typedef GpuPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;
	typedef Link<GpuTimeslice,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(ptrToDeviceSize );
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	PhiloxWrapper rng( site, rngSeed, rngCounter );

	Quaternion<Real> q;
	
	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

		locU.identity();
		
		for( int i=0; i<2; i++ )
			for( int j=i+1; j<3; j++ )
			{
				q[0] = rng.rand()*2.0-1.0;
				q[1] = rng.rand()*2.0-1.0;
				q[2] = rng.rand()*2.0-1.0;
				q[3] = rng.rand()*2.0-1.0;
				
				q.projectSU2();
				locU.rightSubgroupMult( i, j, &q );
			}
 		globUp = locU;
	}
}

}


#endif /* COULOMBKERNELSSU2_HXX_ */
