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
 * 
 * This class contains the kernels and wrappers called
 * by MultiGPU_MPI_Communicator.hxx.
 * 
 * kernels as class members are not supported (even static): 
 * wrap the kernel calls and hide the kernels in namespace.
 * 
 */

#include "../../kernel_launch_bounds.h"
#include "../../GlobalConstants.h"
#include "../../GaugeFixingSubgroupStep.hxx"
#include "../../algorithms/OrUpdate.hxx"
#include "../../algorithms/MicroUpdate.hxx"
#include "../../algorithms/SaUpdate.hxx"
#include "../../algorithms/RandomUpdate.hxx"
#include "../../../lattice/access_pattern/StandardPattern.hxx"
#include "../../../lattice/access_pattern/GpuPatternTimesliceParityPriority.hxx"
#include "../../../lattice/access_pattern/GpuPatternParityPriority.hxx"
#include "../../../lattice/SiteCoord.hxx"
#include "../../../lattice/SiteIndex.hxx"
#include "../../../lattice/Link.hxx"
#include "../../../lattice/SU3.hxx"
#include "../../../lattice/Matrix.hxx"
#include "../../../lattice/rng/PhiloxWrapper.hxx"

#include "./MultiGPU_MPI_LandauKernelsSU3.h"


// kernels:
namespace MPILKSU3
{

template<class Algorithm> inline __device__ void applyOneTimeslice( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, Algorithm algorithm  )
{
	typedef GpuPatternParityPriority< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
// 	SiteIndex<4,FULL_SPLIT> s( DEVICE_CONSTANTS::SIZE_TIMESLICE );
	
	s.nn = nnt;

	const bool updown = threadIdx.x / (NSB*4);
	const short mu = (threadIdx.x % (NSB*4)) / NSB;
	const short id = (threadIdx.x % (NSB*4)) % NSB;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );

	Real* U;
	
	if( updown==1 )
	{
		if( mu!=0 )
		{
			s.setNeighbour(mu,false);
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
	typedef GpuPatternParityPriority< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
// 	SiteIndex<4,FULL_SPLIT> s( DEVICE_CONSTANTS::SIZE_TIMESLICE);
	
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



__global__ void __launch_bounds__(8*NSB,OR_MINBLOCKS) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
{
	OrUpdate overrelax( orParameter );
	applyOneTimeslice( UtUp, UtDw, nnt, parity, overrelax  );
}

__global__ void __launch_bounds__(8*NSB,MS_MINBLOCKS) microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity )
{
	MicroUpdate micro;
	applyOneTimeslice( UtUp, UtDw, nnt, parity, micro );
}

__global__ void __launch_bounds__(8*NSB,SA_MINBLOCKS)  saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	SaUpdate sa( temperature, &rng );
	applyOneTimeslice( UtUp, UtDw, nnt, parity, sa );
}

__global__ void randomTrafo( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, int rngSeed, int rngCounter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, rngSeed, rngCounter );
	RandomUpdate random( &rng );
	applyOneTimeslice( UtUp, UtDw, nnt, parity, random );
}

__global__ void projectSU3( Real* Ut )
{
	typedef GpuPatternParityPriority< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLinkIndex linkUp( Ut, s, mu );
		SU3<TLinkIndex> globUp( linkUp );

		globUp.projectSU3(); // IMPORTANT: Currently this kernel is used for reconstructing third line in the end. Be aware of this when changing something.
	}
}

__global__ void setHot( Real* Ut, int rngSeed, int rngCounter )
{
	typedef GpuPatternParityPriority< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	s.setLatticeIndex( site );
	
	PhiloxWrapper rng( site, rngSeed, rngCounter );

	Quaternion<Real> q;
	
	for( int mu = 0; mu < 4; mu++ )
	{
		TLinkIndex linkUp( Ut, s, mu );
		SU3<TLinkIndex> globUp( linkUp );

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



// wrappers:
// constructor
MultiGPU_MPI_LandauKernelsSU3::MultiGPU_MPI_LandauKernelsSU3()
{
}

// TODO call this function in main or somewhere appropriate
void MultiGPU_MPI_LandauKernelsSU3::initCacheConfig()
{
	cudaFuncSetCacheConfig( MPILKSU3::generateGaugeQualityPerSite, cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( MPILKSU3::randomTrafo, cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( MPILKSU3::orStep, cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( MPILKSU3::microStep, cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( MPILKSU3::saStep, cudaFuncCachePreferL1 );
	cudaFuncSetCacheConfig( MPILKSU3::projectSU3, cudaFuncCachePreferL1 );
}


void MultiGPU_MPI_LandauKernelsSU3::applyOneTimeslice( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, MultiGPU_MPI_AlgorithmOptions algoOptions )
{
	switch( algoOptions.getAlgorithm() )
	{
	case OR:
		MPILKSU3::orStep<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, algoOptions.getOrParameter() );
		break;
	case MS:
		MPILKSU3::microStep<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity );
		break;
	case SA:
		MPILKSU3::saStep<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, algoOptions.getTemperature(), PhiloxWrapper::getNextCounter(), algoOptions.getSeed() );
		break;
	case RT:
		MPILKSU3::randomTrafo<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, PhiloxWrapper::getNextCounter(), algoOptions.getSeed() );
		break;
	default:
		printf("Algorithm type not set to a known value [MultiGPU_MPI_AlgorithmOptions::setAlgorithm(enum AlgoType)]. Exiting\n");
		exit(1);
	}
}

void MultiGPU_MPI_LandauKernelsSU3::projectSU3( int a, int b, cudaStream_t stream, Real* Ut )
{
	MPILKSU3::projectSU3<<<a,b,0,stream>>>( Ut );
}

void MultiGPU_MPI_LandauKernelsSU3::setHot( int a, int b, cudaStream_t stream, Real* Ut, MultiGPU_MPI_AlgorithmOptions algoOptions )
{
	
	MPILKSU3::setHot<<<a,b,0,stream>>>( Ut, PhiloxWrapper::getNextCounter(), algoOptions.getSeed() );
}

void MultiGPU_MPI_LandauKernelsSU3::generateGaugeQualityPerSite( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA )
{
	MPILKSU3::generateGaugeQualityPerSite<<<a,b,0,stream>>>( UtUp, UtDw, nnt, parity, dGff, dA );
}






