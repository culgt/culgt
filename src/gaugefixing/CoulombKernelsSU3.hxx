/*
 * CoulombKernelsSU3.hxx
 *
 *
 *	TODO: simplify the kernels: micro, or, sa are acutally copies of each other with different Algorithms...
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
 */

#ifndef COULOMBKERNELSSU3_HXX_
#define COULOMBKERNELSSU3_HXX_

#include "GlobalConstants.h"
#include "kernel_launch_bounds.h"
#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "../lattice/rng/PhiloxWrapper.hxx"
#include "GaugeFixingSubgroupStep.hxx"
#include "algorithms/SaUpdate.hxx"
#include "algorithms/OrUpdate.hxx"
#include "algorithms/MicroUpdate.hxx"
#include "algorithms/RandomUpdate.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Link.hxx"

// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.
namespace CKSU3
{
static const int Ndim = 4; // TODO why here?
static const int Nc = 3;
__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA );
__global__ void restoreThirdLine( Real* U, lat_index_t* nnt );
__global__ void randomTrafo( Real* UtUp, Real* UtDw,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter );
__global__ void orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter );
__global__ void microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity );
__global__ void saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter );
}

class CoulombKernelsSU3
{
public:
//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );

	// TODO remove static and make the init in constructor
	static void initCacheConfig()
	{
		cudaFuncSetCacheConfig( CKSU3::generateGaugeQualityPerSite, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( CKSU3::restoreThirdLine, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( CKSU3::randomTrafo, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( CKSU3::orStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( CKSU3::microStep, cudaFuncCachePreferL1 );
		cudaFuncSetCacheConfig( CKSU3::saStep, cudaFuncCachePreferL1 );
	}

	static void generateGaugeQualityPerSite( int a, int b, Real *U, double *dGff, double *dA )
	{
		CKSU3::generateGaugeQualityPerSite<<<a,b>>>(U,dGff, dA);
	};
	static double getGaugeQualityPrefactorA()
	{
		return 1./(double)CKSU3::Nc;
	};
	static double getGaugeQualityPrefactorGff()
	{
		return 1./(double)(CKSU3::Nc*(CKSU3::Ndim-1));
	};

	static void restoreThirdLine( int a, int b, Real* U, lat_index_t* nnt )
	{
		CKSU3::restoreThirdLine<<<a,b>>>(U,nnt);
	};
	static void randomTrafo( int a, int b, Real* UtUp, Real* UtDw,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
	{
		CKSU3::randomTrafo<<<a,b>>>( UtUp, UtDw, nnt, parity, rngSeed, rngCounter );
	};
	static void orStep( int a, int b,  Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
	{
		CKSU3::orStep<<<a,b>>>( UtUp, UtDw, nnt, parity, orParameter );
	};
	static void microStep( int a, int b, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity )
	{
		CKSU3::microStep<<<a,b>>>( UtUp, UtDw, nnt, parity );
	};
	static void saStep( int a, int b, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
	{
		CKSU3::saStep<<<a,b>>>( UtUp, UtDw, nnt, parity, temperature, rngSeed, rngCounter);
	};
private:
};





namespace CKSU3
{

__global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA )
{
	typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

	SiteCoord<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE_TIMESLICE);

	int site = blockIdx.x * blockDim.x + threadIdx.x;

	Matrix<Complex<Real>,Nc> locMatSum;
	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);

	Sum.zero();

	// TODO check if there is a faster way to compute DELTA
	for( int mu = 1; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );

		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);

		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		temp.assignWithoutThirdLine( globUp );
//				temp.projectSU3withoutThirdRow();// TODO project here?
//				globUp.assignWithoutThirdLine( temp ); // TODO
		temp.reconstructThirdLine();
		Sum += temp;

		s.setNeighbour(mu,-1);
		TLink linkDw( U, s, mu );
		SU3<TLink> globDw( linkDw );
		temp.assignWithoutThirdLine( globDw );
//				temp.projectSU3withoutThirdRow(); // TODO project here?
//				globDw.assignWithoutThirdLine( temp ); // TODO
		temp.reconstructThirdLine();
		Sum -= temp;
	}

	Sum -= Sum.trace()/Real(3.);

	Matrix<Complex<Real>,Nc> locMatSumHerm;
	SU3<Matrix<Complex<Real>,Nc> > SumHerm(locMatSumHerm);
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	double prec = 0;
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			prec += Sum.get(i,j).abs_squared();
		}
	}

	dA[site] = prec;


	s.setLatticeIndex( site );
	double result = 0;

	Matrix<Complex<Real>,Nc> locTemp;
	SU3<Matrix<Complex<Real>,Nc> > temp(locTemp);
	for( int mu = 1; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );
		temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
		temp.reconstructThirdLine();
		result += temp.trace().x;
	}


//	Matrix<Complex<Real>,Nc> locTemp;
//	SU3<Matrix<Complex<Real>,Nc> > temp(locTemp);
//	Matrix<Complex<double>,Nc > doubleMat;
////	SU3<Matrix<Complex<double>,Nc > > doubleTemp(doubleMat);
//	for( int mu = 1; mu < 4; mu++ )
//	{
//		TLink linkUp( U, s, mu );
//		SU3<TLink> globUp( linkUp );
//		temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
//
////		doubleTemp = temp;
//
//		for( int i = 0; i < 2; i++ )
//			for(int j = 0; j < 3; j++ )
//			{
//				Complex<double> a;
//				a.x = (double)temp.get( i, j ).x;
//				a.y = (double)temp.get( i, j ).y;
//				doubleMat.set( i,j, a);
//			}
//
//
//		Complex<double> a = doubleMat.get(0,1).conj() * doubleMat.get(1,2).conj() - doubleMat.get(0,2).conj() * doubleMat.get(1,1).conj();
//		Complex<double> b = doubleMat.get(0,2).conj() * doubleMat.get(1,0).conj() - doubleMat.get(0,0).conj() * doubleMat.get(1,2).conj();
//		Complex<double> c = doubleMat.get(0,0).conj() * doubleMat.get(1,1).conj() - doubleMat.get(0,1).conj() * doubleMat.get(1,0).conj();
//
//
//		doubleMat.set( 2, 0, a );
//		doubleMat.set( 2, 1, b );
//		doubleMat.set( 2, 2, c );
//
//		result += doubleMat.trace().x;
//	}

	dGff[site] = result;

}

__global__ void restoreThirdLine( Real* U, lat_index_t* nnt )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
	typedef Link<Gpu,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;

//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE_TIMESLICE);
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



template<class Algorithm> inline __device__ void apply( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, Algorithm algorithm  )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink3_2;

	lat_coord_t size[4] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s( size ); // If i give DEVICE_CONSTANTS::SIZE_TIMESLICE instead, register spilling is much higher! Why?
	s.nn = nnt;

	const bool updown = threadIdx.x / (4*NSB);
	const short mu = (threadIdx.x % (4*NSB)) / NSB;
	const short id = (threadIdx.x % (4*NSB)) % NSB;

	int site = blockIdx.x * blockDim.x/8 + id;

	if( parity == true ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );
	if( (mu!=0)&&(updown==true) )
	{
		s.setNeighbour(mu,false);
	}

	Matrix<Complex<Real>,Nc> locMat;
	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

	TLink3_2 link( ((mu==0)&&(updown==true))?(UtDw):(UtUp), s, mu );

	SU3<TLink3_2> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();

	// define the update algorithm
//	OrUpdate overrelax( orParameter );
	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, Algorithm, COULOMB> subgroupStep( &locU, algorithm, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}

__global__ void __launch_bounds__(256,OR_MINBLOCKS) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
{
	OrUpdate overrelax( orParameter );
	apply( UtUp, UtDw, nnt, parity, overrelax );
}



__global__ void __launch_bounds__(256,MS_MINBLOCKS) microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity )
{
	MicroUpdate micro;
	apply( UtUp, UtDw, nnt, parity, micro );
}


__global__ void __launch_bounds__(256,SA_MINBLOCKS) saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	SaUpdate sa( temperature, &rng );
	apply( UtUp, UtDw, nnt, parity, sa );
}

__global__ void randomTrafo( Real* UtUp, Real* UtDw,lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	RandomUpdate random( &rng );
	apply( UtUp, UtDw, nnt, parity, random );
}

}


#endif /* COULOMBKERNELSSU2_HXX_ */
