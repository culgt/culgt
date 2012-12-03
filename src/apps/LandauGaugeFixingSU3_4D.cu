/*
 * test_gaugefixing.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: vogt
 */

#include <iostream>
#include <math.h>
#include <sstream>
#ifndef OSX
#include "malloc.h"
#endif
#include "../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
#include "../lattice/gaugefixing/GaugeFixingStats.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../lattice/gaugefixing/overrelaxation/MicroUpdate.hxx"
//#include "../lattice/gaugefixing/overrelaxation/SrUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/LinkFile.hxx"
//#include "../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/filetype_typedefs.h"
#include "../lattice/gaugefixing/GlobalConstants.hxx"
#include "../lattice/gaugefixing/LandauKernelsSU3.hxx"
#include "../lattice/gaugefixing/CommonKernelsSU3.hxx"
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;



const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;


//typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;


void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s( HOST_CONSTANTS::SIZE);
	s.calculateNeighbourTable( nnt );
}



//// TODO Hack, because cuda5.0 has a bug and the native *= operator does not work
//__device__ Matrix<Complex<Real>,3>& mult( Matrix<Complex<Real>,3>& c, Matrix<Complex<Real>,3>& a, Matrix<Complex<Real>,3>& b )
//{
//	for(int i = 0; i < 3; i++ )
//	{
//		for( int j = 0; j < 3; j++ )
//		{
//			c.mat[i*3+j] = 0;
//			for( int k = 0; k < 3; k++ )
//			{
//				c.mat[i*3+j] += a.mat[i*3+k] * b.mat[k*3+j];
//			}
//		}
//	}
//	return c;
//}
//
//__global__ void calculatePlaquette( Real *U, lat_index_t* nn, double *dPlaquette )
//{
//	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
//	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
//
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteIndex<4,FULL_SPLIT> s(size);
//	s.nn = nn;
//
//	Matrix<Complex<Real>,Nc> matP;
////	SU3<Matrix<Complex<Real>,Nc> > P(matP);
//
//	Matrix<Complex<Real>,Nc> matTemp;
//	SU3<Matrix<Complex<Real>,Nc> > temp(matTemp);
//
//	Matrix<Complex<Real>,Nc> matTemp2;
//
//
//	double localPlaquette = 0;
//
//
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		for( int nu = mu+1; nu < 4; nu++)
//		{
//			for( int i = 0; i < 3; i++ )
//				for( int j =0; j < 3; j++ )
//				{
//					if( i == j ) matP.set( i,j,Complex<Real>(1.0,.0) );
//					else matP.set( i,j,Complex<Real>(.0,.0) );
//				}
//
//			{
//				s.setLatticeIndex( site );
//
//				TLinkIndex link( U, s, mu );
//				SU3<TLinkIndex> globU( link );
//
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//
//				mult(matTemp2, matP, temp.mat );
//				matP = matTemp2;
////				matTemp2 = matP * matTemp;
//			}
//
//			{
//				s.setNeighbour( mu, true );
//
//				TLinkIndex link( U, s, nu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//			{
//				s.setLatticeIndex( site );
//				s.setNeighbour(nu, true );
//
//				TLinkIndex link( U, s, mu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//				temp.hermitian();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//			{
//				s.setLatticeIndex( site );
//
//				TLinkIndex link( U, s, nu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//				temp.hermitian();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//
//
//
//
//			localPlaquette += matP.trace().x;
//		}
//	}
//
//	dPlaquette[site] = localPlaquette/6./3.;
////	dPlaquette[site] = 1;
////	dPlaquette[site] = 0;
//}
//
//__global__ void printPlaquette( double* dPlaquette )
//{
//	const lat_coord_t size[Ndim] = {Nx,Ny,Nz,Nt};
//	SiteCoord<4,FULL_SPLIT> s(size);
//
//	double plaquette = 0;
//	for( int i = 0; i < s.getLatticeSize(); i++ )
//	{
//		plaquette += dPlaquette[i];
//	}
//
//	printf( "\t%E\n", plaquette/double(s.getLatticeSize()) );
//
//}

//Real calculatePolyakovLoopAverage( Real *U )
//{
//	Matrix<Complex<Real>,3> tempMat;
//	SU3<Matrix<Complex<Real>,3> > temp( tempMat );
//	Matrix<Complex<Real>,3> temp2Mat;
//	SU3<Matrix<Complex<Real>,3> > temp2( temp2Mat );
//
//	SiteCoord<Ndim,FULL_SPLIT> s( HOST_CONSTANTS::SIZE );
//
//	Complex<Real> result(0,0);
//
//	for( s[1] = 0; s[1] < s.size[1]; s[1]++ )
//	{
//		for( s[2] = 0; s[2] < s.size[2]; s[2]++ )
//		{
//			for( s[3] = 0; s[3] < s.size[3]; s[3]++ )
//			{
//				temp.identity();
//				temp2.zero();
//
//				for( s[0] = 0; s[0] < s.size[0]; s[0]++ )
//				{
//
//					TLink link( U, s, 0 );
//					SU3<TLink> globU( link );
//
//					temp2 = temp2 + temp*globU;
//
//					temp = temp2;
//					temp2.zero();
//				}
//				result += temp.trace();
//			}
//		}
//	}
//
//	return sqrt(result.x*result.x+result.y*result.y) / (Real)(s.getLatticeSizeTimeslice()*Nc);
//}






int main(int argc, char* argv[])
{
	Chronotimer allTimer;
	allTimer.reset();
	allTimer.start();

	LandauKernelsSU3::initCacheConfig();

	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;

	// Choose device and print device infos
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, options.getDeviceNumber() );
	cudaSetDevice( options.getDeviceNumber() );

	printf("\nDevice %d: \"%s\"\n", options.getDeviceNumber(), deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);


	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


	// TODO maybe we should choose the filetype at compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly;
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt;
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain;


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for configuration
	Real* dU;
	cudaMalloc( &dU, arraySize*sizeof(Real) );

	// host memory for the neighbour table
	lat_index_t* nn = (lat_index_t*)malloc( s.getLatticeSize()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the neighbour table
	lat_index_t *dNn;
	cudaMalloc( &dNn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ) );

	double* dPlaquette;
	cudaMalloc( &dPlaquette, s.getLatticeSize()*sizeof(double) );

	// initialise the neighbour table
	initNeighbourTable( nn );
	// copy neighbour table to device
	cudaMemcpy( dNn, nn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = NSB*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSize()/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call


	GaugeFixingStats<Ndim,Nc,LandauKernelsSU3,AVERAGE> gaugeStats( dU, HOST_CONSTANTS::SIZE );


	// timer to measure kernel times
	Chronotimer kernelTimer;
	kernelTimer.reset();
	kernelTimer.start();

	double orTotalKernelTime = 0; // sum up total kernel time for OR
	long orTotalStepnumber = 0;
	double saTotalKernelTime = 0;

	FileIterator fi( options );

	for( fi.reset(); fi.hasNext(); fi.next() )
	{
		cout << "loading " << fi.getFilename() << " as " << options.getFType() << endl;

		bool loadOk;

		switch( options.getFType() )
		{
		case VOGT:
			lfVogt.reinterpret = options.getReinterpret();
			loadOk = lfVogt.load( s, fi.getFilename(), U );
			break;
		case PLAIN:
			lfPlain.reinterpret = options.getReinterpret();
			loadOk = lfPlain.load( s, fi.getFilename(), U );
			break;
		case HEADERONLY:
			lfHeaderOnly.reinterpret = options.getReinterpret();
			loadOk = lfHeaderOnly.load( s, fi.getFilename(), U );
			break;
		default:
			cout << "Filetype not set to a known value. Exiting";
			exit(1);
		}

		if( !loadOk )
		{
			cout << "Error while loading. Trying next file." << endl;
			break;
		}
		else
		{
			cout << "File loaded." << endl;
		}


		double bestGff = 0.0;
		for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
		{
			// we copy from host in every gaugecopy step to have a cleaner configuration (concerning numerical errors)
			cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

			if( !options.isNoRandomTrafo() ) // I'm an optimist! This should be called isRandomTrafo()!
			{
				LandauKernelsSU3::randomTrafo(numBlocks,threadsPerBlock,dU, dNn, 0, PhiloxWrapper::getNextCounter() );
				LandauKernelsSU3::randomTrafo(numBlocks,threadsPerBlock,dU, dNn, 1, PhiloxWrapper::getNextCounter() );
			}

			// calculate and print the gauge quality
			printf( "i:\t\tgff:\t\tdA:\n");
			gaugeStats.generateGaugeQuality();
			printf( "   \t\t%1.10f\t\t%e\n", gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

			float temperature = options.getSaMax();
			float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getSaSteps(); i++ )
			{
				LandauKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 0, temperature, PhiloxWrapper::getNextCounter() );
				LandauKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 1, temperature, PhiloxWrapper::getNextCounter() );

				for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
				{
					LandauKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 0 );
					LandauKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 1 );
				}


				if( i % options.getCheckPrecision() == 0 )
				{
					CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
	//				projectSU3<<<numBlocks*2,32>>>( dU );

					gaugeStats.generateGaugeQuality();
	//				CudaError::getLastError( "generateGaugeQuality error" );
					printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
				}
				temperature -= tempStep;
			}
			cudaThreadSynchronize();
			kernelTimer.stop();
			saTotalKernelTime += kernelTimer.getTime();

			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{

				LandauKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getOrParameter() );
				LandauKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getOrParameter() );

				if( i % options.getCheckPrecision() == 0 )
				{
					CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
	//				projectSU3<<<numBlocks*2,32>>>( dU );
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
				}

				orTotalStepnumber++;
			}

			cudaThreadSynchronize();
			kernelTimer.stop();
			cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
			orTotalKernelTime += kernelTimer.getTime();




			// reconstruct third line
			CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );

			// check for best copy
			if( gaugeStats.getCurrentGff() > bestGff )
			{
				cout << "FOUND BETTER COPY" << endl;
				bestGff = gaugeStats.getCurrentGff();

				// copy back
				cudaMemcpy( U, dU, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );
			}
			else
			{
				cout << "NO BETTER COPY" << endl;
			}
		}
		
		//saving file
		cout << "saving " << fi.getOutputFilename() << " as " << options.getFType() << endl;

		switch( options.getFType() )
		{
		case VOGT:
			loadOk = lfVogt.save( s, fi.getOutputFilename(), U );
			break;
		case PLAIN:
			loadOk = lfPlain.save( s, fi.getOutputFilename(), U );
			break;
		case HEADERONLY:
			loadOk = lfHeaderOnly.save( s, fi.getOutputFilename(), U );
			break;
		default:
			cout << "Filetype not set to a known value. Exiting";
			exit(1);
		}
	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
//	cout << "total kernel time: " << orTotalKernelTime << " s" << endl;



	long hbFlops = 2176;
	long microFlops = 2118;
	cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


//	long orFlops = 2253;
	long orFlops = 2124; // HV 2012-12-03
	cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

}
