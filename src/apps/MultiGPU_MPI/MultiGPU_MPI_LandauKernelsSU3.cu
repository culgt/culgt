/*
 * CoulombKernelsSU3.hxx
 *
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
 */

#include "../../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
// #include "../../lattice/gaugefixing/GaugeFixingStats.hxx"
#include "../../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../../lattice/access_pattern/GpuLandauPatternParity.hxx"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/SiteIndex.hxx"
#include "../../lattice/Link.hxx"
#include "../../lattice/SU3.hxx"
#include "../../lattice/Matrix.hxx"
#include "../../lattice/gaugefixing/GlobalConstants.hxx"
// #include "../../util/rng/PhiloxWrapper.hxx"

#include "./MultiGPU_MPI_LandauKernelsSU3.h"


// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.

// kernels:
namespace MPILKSU3
{

template<class Algorithm> inline __global__ void __launch_bounds__(256,4) applyOneTimeslice( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, Algorithm algorithm  )
{
	typedef GpuLandauPatternParity< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	s.nn = nnt;

	const bool updown = threadIdx.x / NSB4;
	const short mu = (threadIdx.x % NSB4) / NSB;
	const short id = (threadIdx.x % NSB4) % NSB;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );

	Real* U;
	
	if( updown==1 )
	{
		if( mu!=0 )
		{
			s.setNeighbour(mu,0);
			U=UtUp;
		}
		else
		{
			U=UtDw;
		}
	}
	else
		U=UtUp;

	Matrix<Complex<Real>,Nc> locMat;
	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

	TLinkIndex link( U, s, mu );

	SU3<TLinkIndex> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();


	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, Algorithm, LANDAU> subgroupStep( &locU, algorithm, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}

__global__ void generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA )
{
	typedef GpuLandauPatternParity< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	s.nn = nnt;
	
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	int resid = site;
	if( parity == 1 ) site += Nx*Ny*Nz/2;
// 	if( site >= Nx*Ny*Nz ) return; //important in case Nx^3 is not power of 2

	Matrix<Complex<Real>,Nc> locMatSum;
	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);
	Sum.zero();
	double result = 0;


	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);
		TLinkIndex linkUp( UtUp, s, mu );
		SU3<TLinkIndex> globUp( linkUp );
		
		temp.assignWithoutThirdLine( globUp );
		temp.reconstructThirdLine();
		result += temp.trace().x;
		Sum += temp;
		
		if( mu==0 )
		{
			TLinkIndex linkDw( UtDw, s, mu );
			SU3<TLinkIndex> globDw( linkDw );
			temp.assignWithoutThirdLine( globDw );
			temp.reconstructThirdLine();
		}
		else
		{
			s.setNeighbour(mu,false);
			TLinkIndex linkDw( UtUp, s, mu );
			SU3<TLinkIndex> globDw( linkDw );
			temp.assignWithoutThirdLine( globDw );
			temp.reconstructThirdLine();
		}
		Sum -= temp;
	}
	dGff[resid] = result;

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

	dA[resid] = prec;

}

__global__ void restoreThirdLine( Real* U, lat_index_t* nnt )
{
// 	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
// 	typedef Link<Gpu,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;
// 
// //	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteIndex<Ndim,FULL_SPLIT> s(DEVICE_CONSTANTS::SIZE);
// 	s.nn = nnt;
// 
// 	int site = blockIdx.x * blockDim.x + threadIdx.x;
// 
// 	s.setLatticeIndex( site );
// 
// 	for( int mu = 0; mu < 4; mu++ )
// 	{
// 		TLink link( U, s, mu );
// 		SU3<TLink> glob( link );
// 		glob.projectSU3();
// 	}
}


__global__ void randomTrafo( Real* U,lat_index_t* nnt, bool parity, int counter )
{
	// TODO
}

__global__ void __launch_bounds__(256,4) orStep( Real* U, lat_index_t* nnt, bool parity, float orParameter )
{
// 	OrUpdate overrelax( orParameter );
// 	apply( U, nnt, parity, overrelax );
}

__global__ void __launch_bounds__(256,4) microStep( Real* U, lat_index_t* nnt, bool parity )
{
// 	MicroUpdate micro;
// 	apply( U, nnt, parity, micro );
}

__global__ void saStep( Real* U, lat_index_t* nnt, bool parity, float temperature, int counter )
{
// 	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, 12345, counter );
// 	SaUpdate sa( temperature, &rng );
// 	apply( U, nnt, parity, sa );
}

}


//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );


// wrappers:
// constructor
MultiGPU_MPI_LandauKernelsSU3::MultiGPU_MPI_LandauKernelsSU3()
{
}

// TODO call this function in main or somewhere appropriate
void MultiGPU_MPI_LandauKernelsSU3::initCacheConfig()
{
// 	cudaFuncSetCacheConfig( MPILKSU3::generateGaugeQualityPerSite, cudaFuncCachePreferL1 );
// 	cudaFuncSetCacheConfig( MPILKSU3::restoreThirdLine, cudaFuncCachePreferL1 );
// 	cudaFuncSetCacheConfig( MPILKSU3::randomTrafo, cudaFuncCachePreferL1 );
// 	cudaFuncSetCacheConfig( MPILKSU3::orStep, cudaFuncCachePreferL1 );
// 	cudaFuncSetCacheConfig( MPILKSU3::microStep, cudaFuncCachePreferL1 );
// 	cudaFuncSetCacheConfig( MPILKSU3::saStep, cudaFuncCachePreferL1 );
}


void MultiGPU_MPI_LandauKernelsSU3::applyOneTimeslice( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, enum AlgoType algorithm )
{
	//TODO maybe somewhere else:
// 	cudaFuncSetCacheConfig( MPILKSU3::applyOneTimeslice, cudaFuncCachePreferL1 );
	
	static OrUpdate overrelax(1.35); //TODO parameters!
	
	switch( algorithm )
	{
	case OR:
		MPILKSU3::applyOneTimeslice<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, overrelax  );
		break;
	case SR:
		//TODO SR
		printf("Algorithm type not set to a known value. Exiting\n");
		exit(1);
		break;
	case SA:
		//TODO SA
		printf("Algorithm type not set to a known value. Exiting\n");
		exit(1);
		break;
	default:
		printf("Algorithm type not set to a known value. Exiting\n");
		exit(1);
	}
}

void MultiGPU_MPI_LandauKernelsSU3::generateGaugeQualityPerSite( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA )
{
	MPILKSU3::generateGaugeQualityPerSite<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, dGff, dA );
};
// double MultiGPU_MPI_LandauKernelsSU3::getGaugeQualityPrefactorA()
// {
// 	return 1./(double)MPILKSU3::Nc;
// };
// double MultiGPU_MPI_LandauKernelsSU3::getGaugeQualityPrefactorGff()
// {
// 	return 1./(double)(MPILKSU3::Nc*MPILKSU3::Ndim);
// };
// 
// void MultiGPU_MPI_LandauKernelsSU3::restoreThirdLine( int a, int b, Real* U, lat_index_t* nnt )
// {
// 	MPILKSU3::restoreThirdLine<<<a,b>>>(U,nnt);
// };
// void MultiGPU_MPI_LandauKernelsSU3::randomTrafo( int a, int b, Real* U,lat_index_t* nnt, bool parity, int counter )
// {
// 	MPILKSU3::randomTrafo<<<a,b>>>( U, nnt, parity, counter );
// };
// void MultiGPU_MPI_LandauKernelsSU3::orStep( int a, int b,  Real* U, lat_index_t* nnt, bool parity, float orParameter )
// {
// 	MPILKSU3::orStep<<<a,b>>>( U, nnt, parity, orParameter );
// };
// void MultiGPU_MPI_LandauKernelsSU3::microStep( int a, int b, Real* U, lat_index_t* nnt, bool parity )
// {
// 	MPILKSU3::microStep<<<a,b>>>( U, nnt, parity );
// };
// void MultiGPU_MPI_LandauKernelsSU3::saStep( int a, int b, Real* U, lat_index_t* nnt, bool parity, float temperature, int counter )
// {
// 	MPILKSU3::saStep<<<a,b>>>( U, nnt, parity, temperature, counter);
// };




