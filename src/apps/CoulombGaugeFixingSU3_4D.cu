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
#include "../lattice/gaugefixing/simulated_annealing/SaUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/LinkFile.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
#include "../lattice/filetypes/filetype_typedefs.h"
#include "../util/datatype/lattice_typedefs.h"
#include "../util/rng/PhiloxWrapper.hxx"
//#include <boost/program_options/parsers.hpp>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/options_description.hpp>
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"
#include "../lattice/gaugefixing/GlobalConstants.hxx"
#include "../util/cuda/CudaError.hxx"
#include "../lattice/gaugefixing/CoulombKernelsSU3.hxx"
#include "../lattice/gaugefixing/CommonKernelsSU3.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;



// lattice setup
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef GpuCoulombPattern<SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;

//typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;
//typedef Link<GpuTimeslice,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink3;

void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE_TIMESLICE);
	s.calculateNeighbourTable( nnt );
}

//__global__ void projectSU3DP( Real* U )
//{
//	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
//	SiteCoord<4,FULL_SPLIT> s(size);
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	s.setLatticeIndex( site );
//
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		TLink3 linkUp( U, s, mu );
//		SU3<TLink3> globUp( linkUp );
//
//		Matrix<complex,Nc> locMat;
//		SU3<Matrix<complex,Nc> > temp(locMat);
//
//		temp.assignWithoutThirdLine( globUp );
//
//		Matrix<Complex<double>,Nc > doubleMat;
//
//		for( int i = 0; i < 2; i++ )
//			for(int j = 0; j < 3; j++ )
//			{
//				Complex<double> a;
//				a.x = (double)temp.get( i, j ).x;
//				a.y = (double)temp.get( i, j ).y;
//				doubleMat.set( i,j, a);
//			}
//		double abs_u = 0, abs_v = 0;
//		Complex<double> sp(0.,0.);
//
//		// normalize first row
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			abs_u += doubleMat.get( 0, i ).abs_squared();
//		}
//
//		abs_u = sqrt(abs_u);
//
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			doubleMat.set( 0, i, doubleMat.get(0,i)/abs_u );
//		}
//
//		// orthogonalize second row
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			sp += doubleMat.get( 1,i ) * doubleMat.get( 0, i ).conj();
//		}
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			doubleMat.set( 1, i, doubleMat.get(1,i) - doubleMat.get( 0, i)*sp );
//		}
//
//		// normalize second row
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			abs_v += doubleMat.get( 1, i ).abs_squared();
//		}
//		abs_v = sqrt(abs_v);
//		for( lat_group_dim_t i = 0; i < 3; i++ )
//		{
//			doubleMat.set( 1, i, doubleMat.get(1,i)/abs_v );
//		}
//
//
//		for( int i = 0; i < 2; i++ )
//			for(int j = 0; j < 3; j++ )
//			{
//				Complex<Real> a;
//				a.x = (Real)doubleMat.get( i, j ).x;
//				a.y = (Real)doubleMat.get( i, j ).y;
//				temp.set( i,j, a);
//			}
//
//		globUp.assignWithoutThirdLine( temp );
//	}
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

	CoulombKernelsSU3::initCacheConfig();

	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;

	// Choose device and print device infos
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	printf("\nDevice %d: \"%s\"\n", 0, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);


	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);

	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly;
	lfHeaderOnly.reinterpret = options.getReinterpret();
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt;
	lfVogt.reinterpret = options.getReinterpret();
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain;
	lfPlain.reinterpret = options.getReinterpret();


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for timeslice t
	Real* dUtUp;
	cudaMalloc( &dUtUp, timesliceArraySize*sizeof(Real) );

	// device memory for timeslice t-1
	Real* dUtDw;
	cudaMalloc( &dUtDw, timesliceArraySize*sizeof(Real) );

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt;
	cudaMalloc( &dNnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof( lat_index_t ) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	// copy neighbour table to device
	cudaMemcpy( dNnt, nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = 32*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSizeTimeslice()/2/32; // // half of the lattice sites (a parity) are updated in a kernel call


	GaugeFixingStats<Ndim,Nc,CoulombKernelsSU3,AVERAGE> gaugeStats( dUtUp, HOST_CONSTANTS::SIZE_TIMESLICE );

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
		cout << "loading " << fi.getFilename() << " as " << options.getFType() << endl; // TODO why is FileType printed as a number (not as text according to filetype_typedefs.h)

		bool loadOk;

		switch(  options.getFType() )
		{
		case VOGT:
			loadOk = lfVogt.load( s, fi.getFilename(), U );
			break;
		case PLAIN:
			loadOk = lfPlain.load( s, fi.getFilename(), U );
			break;
		case HEADERONLY:
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

		for( int t = 0; t < s.size[0]; t++ )
		{
			int tDw = (t > 0)?(t-1):(s.size[0]-1); // calculating t-1 (periodic boundaries)

			double bestGff = 0.0;
			for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
			{
				// copying timeslice t ...
				cudaMemcpy( dUtUp, &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
				// ... and t-1 to device
				cudaMemcpy( dUtDw, &U[tDw*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
				// TODO it is not necessary to copy the (t-1) again for t>0, simply swap pointers on device side...

				if( !options.isNoRandomTrafo() ) // I'm an optimist! This should be called isRandomTrafo()!
				{
					CoulombKernelsSU3::randomTrafo(numBlocks,threadsPerBlock, dUtUp, dUtDw, dNnt, 0, PhiloxWrapper::getNextCounter() );
					CoulombKernelsSU3::randomTrafo(numBlocks,threadsPerBlock, dUtUp, dUtDw, dNnt, 1, PhiloxWrapper::getNextCounter() );
				}

				// calculate and print the gauge quality
				printf( "i:\t\tgff:\t\tdA:\n");
				gaugeStats.generateGaugeQuality();


				float temperature = options.getSaMax();
				float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

				kernelTimer.reset();
				kernelTimer.start();
				for( int i = 0; i < options.getSaSteps(); i++ )
				{


					CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, temperature, PhiloxWrapper::getNextCounter() );
					CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, temperature, PhiloxWrapper::getNextCounter() );

					for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
					{
						CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0 );
						CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1 );
					}


					if( i % options.getCheckPrecision() == 0 )
					{
						CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
						CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );

						gaugeStats.generateGaugeQuality();
						CudaError::getLastError( "generateGaugeQuality error" );
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

					CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, options.getOrParameter() );
					CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, options.getOrParameter() );

					if( i % options.getCheckPrecision() == 0 )
					{
						CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
						CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );

						gaugeStats.generateGaugeQuality();
						printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

						if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
					}

					orTotalStepnumber++;
				}

				cudaThreadSynchronize();
				kernelTimer.stop();
				orTotalKernelTime += kernelTimer.getTime();

				// reconstruct third line before copy back
				CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtUp, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );
				CommonKernelsSU3::projectSU3( numBlocks*2, 32, dUtDw, HOST_CONSTANTS::getPtrToDeviceSizeTimeslice() );

				if( gaugeStats.getCurrentGff() > bestGff )
				{
					cout << "FOUND BETTER COPY" << endl;
					bestGff = gaugeStats.getCurrentGff();

					// copy back
					cudaMemcpy( &U[t*timesliceArraySize], dUtUp, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
					cudaMemcpy( &U[tDw*timesliceArraySize], dUtDw, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
				}
				else
				{
					cout << "NO BETTER COPY" << endl;
				}
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


	long hbFlops = 2176-18;
	long microFlops = 2118-18;
	cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


//	long orFlops = 2253;
	long orFlops = 2124-18; // HV 2012-12-03
	cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9/(double)s.size[0] << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9/(double)s.size[0] << "GB/s memory throughput." << endl;
}
