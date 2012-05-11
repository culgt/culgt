/*
 * test_gaugefixing.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: vogt
 */

#include <iostream>
#include <math.h>
#include <sstream>
#include <malloc.h>
#include "../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/LinkFile.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/HeaderVogt.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;

#ifdef _X_
const lat_coord_t Nx = _X_;
#else
#error "Define X (the lattice size in x-direction)"
#endif
#ifdef _Y_
const lat_coord_t Ny = _Y_;
#else
const lat_coord_t Ny = _X_;
bool warnY = true; // TODO print the warning
#endif
#ifdef _Z_
const lat_coord_t Nz = _Z_;
#else
const lat_coord_t Nz = _X_;
bool warnZ = true;
#endif
#ifdef _T_
const lat_coord_t Nt = _T_;
#else
#error "Define T (the lattice size in t-direction)"
#endif

const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,false>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;


typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;


void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.calculateNeighbourTable( nnt );
}

__global__ void printGaugeQuality( Real* dGff )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteCoord<4,true> s(size);

	Real gff = 0;
	for( int i = 0; i < s.getLatticeSize(); i++ )
	{
		gff+= dGff[i];
	}

	printf( "gff: %1.10f\n", gff/Real(s.getLatticeSize())/3./4. );

}

__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteCoord<4,true> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		globUp.projectSU3();
	}
}

__global__ void generateGaugeQuality( Real *U, Real *dGff )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteCoord<4,true> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	// TODO calculate DELTA

//	for( int mu = 1; mu < 4; mu++ )
//	{
//		TLink3 linkUp( U, s, mu ); // TODO s should be passed by reference
//		SU3<TLink3> globUp( linkUp );
//
//		s.setNeighbour(mu,-1); // TODO but then this is not possible
//		TLink3 linkDw( U, s, mu );
//		SU3<TLink3> globDw( linkDw );
//
//		Matrix<complex,Nc> locMatUp;
//		SU3<Matrix<complex,Nc> > Aup(locMatUp);
//
//		Matrix<complex,Nc> locMatDw;
//		SU3<Matrix<complex,Nc> > Adw(locMatDw);
//
//		Aup = globUp - globUp.hermitian();
//		Aup = globDw - globDw.hermitian();
//
//		Aup /= complex(0,2);
//		Adw /= complex(0,2);
//
//		complex trUp = Aup.trace();
//		complex trDw = ADw.trace();
//
//		Aup -= complex(0,Real(1/3))*trUp;
//		Adw -= complex(0,Real(1/3))*trDw;
//	}



	Real result = 0;

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );
		result += globUp.trace().x;
	}

	dGff[site] = result;
}



__global__ void __launch_bounds__(256,4) orStep( Real* U, lat_index_t* nn, bool parity, float orParameter )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.nn = nn;

	const bool updown = threadIdx.x / 128;
	const short mu = (threadIdx.x % 128) / 32;
	const short id = (threadIdx.x % 128) % 32;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );
	if( updown==1 )
	{
		s.setNeighbour(mu,0);
	}

//	if(id == 0) printf("bin in or\n");

	Matrix<complex,Nc> locMat;
	SU3<Matrix<complex,Nc> > locU(locMat);

	TLinkIndex link( U, s, mu );

	SU3<TLinkIndex> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();

	// define the update algorithm
	OrUpdate overrelax( orParameter );
	GaugeFixingSubgroupStep<SU3<Matrix<complex,Nc> >, OrUpdate, LANDAU> subgroupStep( &locU, overrelax, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<complex,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}




Real calculatePolyakovLoopAverage( Real *U )
{
	Matrix<complex,3> tempMat;
	SU3<Matrix<complex,3> > temp( tempMat );
	Matrix<complex,3> temp2Mat;
	SU3<Matrix<complex,3> > temp2( temp2Mat );

	SiteCoord<Ndim,true> s( size );

	complex result(0,0);

	for( s[1] = 0; s[1] < s.size[1]; s[1]++ )
	{
		for( s[2] = 0; s[2] < s.size[2]; s[2]++ )
		{
			for( s[3] = 0; s[3] < s.size[3]; s[3]++ )
			{
				temp.identity();
				temp2.zero();

				for( s[0] = 0; s[0] < s.size[0]; s[0]++ )
				{

					TLink link( U, s, 0 );
					SU3<TLink> globU( link );

					temp2 = temp2 + temp*globU;

					temp = temp2;
					temp2.zero();
				}
				result += temp.trace();
			}
		}
	}

	return sqrt(result.x*result.x+result.y*result.y) / (Real)(s.getLatticeSizeTimeslice()*Nc);
}






int main(int argc, char* argv[])
{

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	printf("\nDevice %d: \"%s\"\n", 0, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);

	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,true> s(size);
	LinkFile<HeaderVogt, void, Standard, Gpu, SiteCoord<4,true> > lf;


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for configuration
	Real* dU;
	cudaMalloc( &dU, arraySize*sizeof(Real) );

	// device memory for collecting the parts of the gauge fixing functional
	Real *dGff;
	cudaMalloc( &dGff, s.getLatticeSize()*sizeof(Real) );

	// host memory for the neighbour table
	lat_index_t* nn = (lat_index_t*)malloc( s.getLatticeSize()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNn;
	cudaMalloc( &dNn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nn );
	// copy neighbour table to device
	cudaMemcpy( dNn, nn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = 32*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSize()/2/32; // // half of the lattice sites (a parity) are updated in a kernel call

	allTimer.start();

	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );


	for( int i = 0; i < 1; i++ )
	{

		stringstream filename(stringstream::out);
		filename << "/home/vogt/configs/STUDIENARBEIT/N32/config_n32t32beta570_sp" << setw( 4 ) << setfill( '0' ) << i << ".vogt";

		bool loadOk = lf.load( s, filename.str(), U );

		if( !loadOk )
		{
			cout << "Error while loading. Trying next file." << endl;
			break;
		}
		else
		{
			cout << "File loaded." << endl;
		}

		// copying configuration t ...
		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

		// calculate and print the gauge quality
		generateGaugeQuality<<<numBlocks*2,32>>>(dU, dGff );
		printGaugeQuality<<<1,1>>>(dGff);


		float orParameter = 1.7;

		for( int i = 0; i < 5000; i++ )
		{
			orStep<<<numBlocks,threadsPerBlock>>>(dU, dNn, 0, orParameter );
			orStep<<<numBlocks,threadsPerBlock>>>(dU, dNn, 1, orParameter );

			if( i % 100 == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );
				generateGaugeQuality<<<numBlocks*2,32>>>(dU, dGff );
				printGaugeQuality<<<1,1>>>(dGff);
			}
		}
		cudaMemcpy( U, dU, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );

	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << endl;
}
