// /*
//  * MultiGPU_MPI_LandauGaugeFixingSU3_4D.h
//  *
//  *  Created on: Nov. 16, 2012
//  *      Author: vogt&schroeck
//  * 
//  * TODO this file is temporary.
//  */
// 
// #include <iostream>
// #include <math.h>
// #include <sstream>
// #ifndef OSX
// #include "malloc.h"
// #endif
// #include "../../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
// // #include "../../lattice/gaugefixing/GaugeFixingStats.hxx"
// #include "../../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
// #include "../../lattice/access_pattern/StandardPattern.hxx"
// #include "../../lattice/access_pattern/GpuCoulombPattern.hxx"
// #include "../../lattice/access_pattern/GpuLandauPatternParity.hxx"
// #include "../../lattice/SiteCoord.hxx"
// #include "../../lattice/SiteIndex.hxx"
// #include "../../lattice/Link.hxx"
// #include "../../lattice/SU3.hxx"
// #include "../../lattice/Matrix.hxx"
// #include "../../lattice/gaugefixing/GlobalConstants.hxx"
// #include "../../util/rng/PhiloxWrapper.hxx"
// #include "MultiGPU_MPI_LandauGaugeFixingSU3_4D.h"
// 
// using namespace std;
// 
// const lat_dim_t Ndim = 4;
// const short Nc = 3;
// 
// // lattice setup
// const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
// const lat_coord_t sizeTimeslice[Ndim] = {1,Nx,Ny,Nz};
// __constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
// 
// const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
// const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 
// 
// __global__ void __launch_bounds__(8*NSB,128/NSB) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter  )
// {
// 	typedef GpuLandauPatternParity< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
// 	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
// 
// 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteIndex<4,FULL_SPLIT> s(size);
// 	s.nn = nnt;
// 
// 	const bool updown = threadIdx.x / NSB4;
// 	const short mu = (threadIdx.x % NSB4) / NSB;
// 	const short id = (threadIdx.x % NSB4) % NSB;
// 
// 	int site = blockIdx.x * blockDim.x/8 + id;
// 	if( parity == 1 ) site += s.getLatticeSize()/2;
// 
// 	s.setLatticeIndex( site );
// 
// 	Real* U;
// 	
// 	if( updown==1 )
// 	{
// 		if( mu!=0 )
// 		{
// 			s.setNeighbour(mu,0);
// 			U=UtUp;
// 		}
// 		else
// 		{
// 			U=UtDw;
// 		}
// 	}
// 	else
// 		U=UtUp;
// 
// 	Matrix<Complex<Real>,Nc> locMat;
// 	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
// 	TLinkIndex link( U, s, mu );
// 	SU3<TLinkIndex> globU( link );
// 	
// 	// make link local
// 	locU.assignWithoutThirdLine(globU);
// 	locU.reconstructThirdLine();
// 
// 	// define the update algorithm
// 	OrUpdate overrelax( orParameter );
// 	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, OrUpdate, LANDAU> subgroupStep( &locU, overrelax, id, mu, updown );
// 
// 	// do the subgroup iteration
// 	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );
// 
// 	// copy link back
// 	globU.assignWithoutThirdLine(locU);
// 	//globU=locU; //TODO with or without 3rd line?
// }
// 
// 
// __global__ void generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, double *dGff, double *dA )
// {
// 	typedef GpuLandauPatternParity< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
// 	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
// 
// 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteIndex<4,FULL_SPLIT> s(size);
// 	s.nn = nnt;
// 	
// 	int site = blockIdx.x * blockDim.x + threadIdx.x;
// 	if( site >= Nx*Ny*Nz ) return; //important in case Nx^3 is not power of 2
// 
// 	Matrix<Complex<Real>,Nc> locMatSum;
// 	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);
// 	Sum.zero();
// 	double result = 0;
// 
// 
// 	for( int mu = 0; mu < 4; mu++ )
// 	{
// 		s.setLatticeIndex( site );
// 		Matrix<Complex<Real>,Nc> locMat;
// 		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);
// 		TLinkIndex linkUp( UtUp, s, mu );
// 		SU3<TLinkIndex> globUp( linkUp );
// 		
// 		temp.assignWithoutThirdLine( globUp );
// 		temp.reconstructThirdLine();
// 		result += temp.trace().x;
// 		Sum += temp;
// 		
// 		if( mu==0 )
// 		{
// 			TLinkIndex linkDw( UtDw, s, mu );
// 			SU3<TLinkIndex> globDw( linkDw );
// 			temp.assignWithoutThirdLine( globDw );
// 			temp.reconstructThirdLine();
// 		}
// 		else
// 		{
// 			s.setNeighbour(mu,false);
// 			TLinkIndex linkDw( UtUp, s, mu );
// 			SU3<TLinkIndex> globDw( linkDw );
// 			temp.assignWithoutThirdLine( globDw );
// 			temp.reconstructThirdLine();
// 		}
// 		Sum -= temp;
// 	}
// 	dGff[site] = result;
// 
// 	Sum -= Sum.trace()/Real(3.);
// 
// 	Matrix<Complex<Real>,Nc> locMatSumHerm;
// 	SU3<Matrix<Complex<Real>,Nc> > SumHerm(locMatSumHerm);
// 	SumHerm = Sum;
// 	SumHerm.hermitian();
// 
// 	Sum -= SumHerm;
// 
// 	double prec = 0;
// 	for( int i = 0; i < 3; i++ )
// 	{
// 		for( int j = 0; j < 3; j++ )
// 		{
// 			prec += Sum.get(i,j).abs_squared();
// 		}
// 	}
// 
// 	dA[site] = prec;
// }
// 
// 
// __global__ void averageGaugeQuality( double* dGff, double* dA )
// {
// 	double gff = 0;
// 	double A = 0;
// 	for( int i = 0; i < Nx*Ny*Nz; i++ )
// 	{
// 		gff+= dGff[i];
// 		A += dA[i];
// 	}
// 	
// 	dGff[0] = gff/double(Nx*Ny*Nz)/4./3.;
// 	dA[0]   = A/double(Nx*Ny*Nz)/3.;
// }
// 
// 
// __global__ void set_hot( Real* U, int counter)
// {
// 	typedef GpuLandauPatternParity< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
// 	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
// 
// 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteIndex<4,FULL_SPLIT> s(size);
// 	int site = blockIdx.x * blockDim.x + threadIdx.x;
// 	s.setLatticeIndex( site );
// 	
// 	PhiloxWrapper rng( site, 123, counter );
// 
// 	Quaternion<Real> q;
// 	
// 	for( int mu = 0; mu < 4; mu++ )
// 	{
// 		TLinkIndex linkUp( U, s, mu );
// 		SU3<TLinkIndex> globUp( linkUp );
// 
// 		Matrix<Complex<Real>,Nc> locMat;
// 		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
// 
// 		locU.identity();
// 		
// 		for( int i=0; i<2; i++ )
// 			for( int j=i+1; j<3; j++ )
// 			{
// 				q[0] = rng.rand()*2.0-1.0;
// 				q[1] = rng.rand()*2.0-1.0;
// 				q[2] = rng.rand()*2.0-1.0;
// 				q[3] = rng.rand()*2.0-1.0;
// 				
// 				q.projectSU2();
// 				locU.rightSubgroupMult( i, j, &q );
// 			}
//  		globUp = locU;
// 	}
// }
// 
// // __global__ void projectSU3( Real* U )
// // {
// // 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// // 	SiteCoord<4,FULL_SPLIT> s(size);
// // 	int site = blockIdx.x * blockDim.x + threadIdx.x;
// // 
// // 	s.setLatticeIndex( site );
// // 
// // 	for( int mu = 0; mu < 4; mu++ )
// // 	{
// // 		TLink linkUp( U, s, mu );
// // 		SU3<TLink> globUp( linkUp );
// // 
// // 
// // 		Matrix<Complex<Real>,Nc> locMat;
// // 		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
// // 		
// // 		locU.assignWithoutThirdLine(globUp);
// // 		locU.projectSU3withoutThirdRow();
// // 		globUp.assignWithoutThirdLine(locU);
// // 	}
// // }
// 
// 
// //---------------- kernel wrapper functions etc. ---------------------
// 
// void initDevice( const int device )
// {
// 	cudaSetDevice(device);
// 	cudaDeviceProp deviceProp;
// 	cudaGetDeviceProperties(&deviceProp, device);
// 	printf("\nDevice %d: \"%s\"\n", device, deviceProp.name);
// 	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);
// }
// 
// void _set_hot( Real* U, int counter)
// {
// 	int threadsPerBlock = 32*8;
// 	int numBlocks = Nx*Ny*Nz/2/32;
// 
// 	set_hot<<<numBlocks*2,32>>>( U, counter ); 
// }
// 
// void _generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, double *dGff, double *dA )
// {
// 	generateGaugeQualityPerSite<<<Nx*Nx*Nx/32,32>>>( UtUp, UtDw, nnt, dGff, dA );
// }
// 
// void _averageGaugeQuality( double* dGff, double* dA )
// {
// 	averageGaugeQuality<<<1,1>>>( dGff, dA );
// }
// 
// void _orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter, cudaStream_t stream  )
// {
// 	int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
// 	int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call
// 
// 	orStep<<<numBlocks,threadsPerBlock,0,stream>>>( UtUp, UtDw, nnt, parity, orParameter );
// }
// 
// //--------------------------------------------------------------
