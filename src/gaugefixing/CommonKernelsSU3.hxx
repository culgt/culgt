/*
 * CommonKernelsSU3.hxx
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
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

private:
};





namespace COMKSU3
{

__global__ void projectSU3( Real *U, lat_coord_t* ptrToDeviceSize )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;
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

}


#endif /* COULOMBKERNELSSU2_HXX_ */
