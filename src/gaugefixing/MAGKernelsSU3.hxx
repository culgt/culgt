/*
 * CoulombKernelsSU3.hxx
 *
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
 */

#ifndef MAGKERNELSSU3_HXX_
#define MAGKERNELSSU3_HXX_

#include "GlobalConstants.h"
#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "../lattice/rng/PhiloxWrapper.hxx"
#include "GaugeFixingSubgroupStep.hxx"
#include "algorithms/SaUpdate.hxx"
#include "algorithms/SrUpdate.hxx"
#include "algorithms/OrUpdate.hxx"
#include "algorithms/MicroUpdate.hxx"
#include "algorithms/RandomUpdate.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Link.hxx"


// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.
namespace MAGKSU3
{
static const int Ndim = 4;
static const int Nc = 3;
__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA );
__global__ void restoreThirdLine( Real* U, lat_index_t* nnt );
__global__ void randomTrafo( Real* U,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter );
__global__ void orStep( Real* U, lat_index_t* nnt, bool parity, float orParameter );
__global__ void srStep( Real* U, lat_index_t* nnt, bool parity, float srParameter, int rngSeed, int rngCounter );
__global__ void microStep( Real* U, lat_index_t* nnt, bool parity );
__global__ void saStep( Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter );
}

class MAGKernelsSU3
{
public:
//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );

	// TODO remove static and make the init in constructor
	static void initCacheConfig()
	{
		cudaFuncSetCacheConfig( MAGKSU3::generateGaugeQualityPerSite, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::restoreThirdLine, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::randomTrafo, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::orStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::srStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::microStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( MAGKSU3::saStep, cudaFuncCachePreferL1 );
	}

	static void generateGaugeQualityPerSite( int a, int b, Real *U, double *dGff, double *dA )
	{
		MAGKSU3::generateGaugeQualityPerSite<<<a,b>>>(U,dGff, dA);
	};
	static double getGaugeQualityPrefactorA()
	{
		return 1./(double)MAGKSU3::Nc;
	};
	static double getGaugeQualityPrefactorGff()
	{
		return 1./(double)(MAGKSU3::Nc*MAGKSU3::Ndim);
	};

	static void restoreThirdLine( int a, int b, Real* U, lat_index_t* nnt )
	{
		MAGKSU3::restoreThirdLine<<<a,b>>>(U,nnt);
	};
	static void randomTrafo( int a, int b, Real* U,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
	{
		MAGKSU3::randomTrafo<<<a,b>>>( U, nnt, parity, rngSeed, rngCounter );
	};
	static void orStep( int a, int b,  Real* U, lat_index_t* nnt, bool parity, float orParameter )
	{
		MAGKSU3::orStep<<<a,b>>>( U, nnt, parity, orParameter );
	};
	static void srStep( int a, int b,  Real* U, lat_index_t* nnt, bool parity, float srParameter, int rngSeed, int rngCounter )
	{
		MAGKSU3::srStep<<<a,b>>>( U, nnt, parity, srParameter, rngSeed, rngCounter );
	};
	static void microStep( int a, int b, Real* U, lat_index_t* nnt, bool parity )
	{
		MAGKSU3::microStep<<<a,b>>>( U, nnt, parity );
	};
	static void saStep( int a, int b, Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
	{
		MAGKSU3::saStep<<<a,b>>>( U, nnt, parity, temperature, rngSeed, rngCounter );
	};
private:
};





namespace MAGKSU3
{

__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA )
{
	typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

	SiteCoord<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	Quaternion<Real> u,v;
	Complex<double> X[3];

	dGff[site]=0.0;

	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		s.setNeighbour(mu,-1);
		TLink linkDw( U, s, mu );
		SU3<TLink> globDw( linkDw );

		// action:
		for( int j=0; j<3; j++ )
			dGff[site] += linkUp.get(j,j).abs_squared();

//		// precision:
//		for( int i=0; i<2; i++ )
//			for( int j=i+1; j<3; j++ )
//			{
//				u = globUp.getSubgroupQuaternion( i, j );
//				u.projectSU2();
//
//				v = globDw.getSubgroupQuaternion( i, j );
//				v.projectSU2();
//
//				X[i+j-1].x += -u[0]*u[2] + u[1]*u[3] + v[0]*v[2] + v[1]*v[3];
//				X[i+j-1].y += -u[0]*u[1] - u[2]*u[3] + v[0]*v[1] - v[2]*v[3];
//			}
		for( int i=0; i<2; i++ )
			for( int j=i+1; j<3; j++ )
			{
				u = globUp.getSubgroupQuaternion( i, j );
//				u.projectSU2();

				v = globDw.getSubgroupQuaternion( i, j );
//				v.projectSU2();

				X[i+j-1].x += -u[0]*u[2] + u[1]*u[3] + v[0]*v[2] + v[1]*v[3];
				X[i+j-1].y += -u[0]*u[1] - u[2]*u[3] + v[0]*v[1] - v[2]*v[3];
			}

	}
	dA[site] = X[0].abs_squared()+X[1].abs_squared()+X[2].abs_squared();







//	typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
//	typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;
//
//	SiteCoord<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	Matrix<complex,Nc> locMatSum;
//	SU3<Matrix<complex,Nc> > Sum(locMatSum);
//
//	Sum.zero();
//
//	// TODO check if there is a faster way to compute DELTA
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		s.setLatticeIndex( site );
//
//		Matrix<complex,Nc> locMat;
//		SU3<Matrix<complex,Nc> > temp(locMat);
//
//		TLink linkUp( U, s, mu );
//		SU3<TLink> globUp( linkUp );
//
//		temp.assignWithoutThirdLine( globUp );
////				temp.projectSU3withoutThirdRow();// TODO project here?
////				globUp.assignWithoutThirdLine( temp ); // TODO
//		temp.reconstructThirdLine();
//		Sum += temp;
//
//		s.setNeighbour(mu,-1);
//		TLink linkDw( U, s, mu );
//		SU3<TLink> globDw( linkDw );
//		temp.assignWithoutThirdLine( globDw );
////				temp.projectSU3withoutThirdRow(); // TODO project here?
////				globDw.assignWithoutThirdLine( temp ); // TODO
//		temp.reconstructThirdLine();
//		Sum -= temp;
//	}
//
//	Sum -= Sum.trace()/Real(3.);
//
//	Matrix<complex,Nc> locMatSumHerm;
//	SU3<Matrix<complex,Nc> > SumHerm(locMatSumHerm);
//	SumHerm = Sum;
//	SumHerm.hermitian();
//
//	Sum -= SumHerm;
//
//	double prec = 0;
//
//	for( int i = 0; i < 3; i++ )
//	{
//		for( int j = 0; j < 3; j++ )
//		{
//			prec += Sum.get(i,j).abs_squared();
//		}
//	}
//	dA[site] = prec;
//
//
//	s.setLatticeIndex( site );
//	double result = 0;
//
//	Matrix<complex,Nc> locTemp;
//	SU3<Matrix<complex,Nc> > temp(locTemp);
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		TLink linkUp( U, s, mu );
//		SU3<TLink> globUp( linkUp );
//		temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
//		temp.reconstructThirdLine();
//		result += temp.trace().x;
//	}
//
//	dGff[site] = result;

}

__global__ void restoreThirdLine( Real* U, lat_index_t* nnt )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
	s.nn = nnt;

	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink link( U, s, mu );
		SU3<TLink> glob( link );
		glob.projectSU3();
	}
}



template<class Algorithm> inline __device__ void apply( Real* U, lat_index_t* nn, bool parity, Algorithm algorithm  )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	s.nn = nn;

	const bool updown = threadIdx.x / NSB4;
	const short mu = (threadIdx.x % NSB4) / NSB;
	const short id = (threadIdx.x % NSB4) % NSB;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );
	if( updown==1 )
	{
		s.setNeighbour(mu,false);
	}

//	if(id == 0) printf("bin in or\n");

	Matrix<Complex<Real>,Nc> locMat;
	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

	TLinkIndex link( U, s, mu );

	SU3<TLinkIndex> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();


	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, Algorithm, MAG> subgroupStep( &locU, algorithm, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}

__global__ void __launch_bounds__(256,4) orStep( Real* U, lat_index_t* nnt, bool parity, float orParameter )
{
	OrUpdate overrelax( orParameter );
	apply( U, nnt, parity, overrelax );
}

__global__ void __launch_bounds__(256,4) microStep( Real* U, lat_index_t* nnt, bool parity )
{
	MicroUpdate micro;
	apply( U, nnt, parity, micro );
}

__global__ void saStep( Real* U, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	SaUpdate sa( temperature, &rng );
	apply( U, nnt, parity, sa );
}

__global__ void __launch_bounds__(256,4) srStep( Real* U, lat_index_t* nnt, bool parity, float srParameter, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	SrUpdate sr( srParameter, &rng );
	apply( U, nnt, parity, sr );
}

__global__ void __launch_bounds__(256,4) randomTrafo( Real* U,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	RandomUpdate random( &rng );
	apply( U, nnt, parity, random );
}

}


#endif /* COULOMBKERNELSSU2_HXX_ */
