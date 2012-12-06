/*
 * test_gaugefixing.cpp
 *
 *  Created on: Sep. 10, 2012
 *      Author: vogt & schroeck
 */

#include <iostream>
#include <math.h>
#include <sstream>
#include <malloc.h>
#include "../GlobalConstants.h"
#include "../GaugeFixingStats.hxx"
#include "../MAGKernelsSU3.hxx"
#include "../CommonKernelsSU3.hxx"
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/SiteIndex.hxx"
#include "../../lattice/LinkFile.hxx"
#include "../../util/timer/Chronotimer.h"
#include "../../lattice/filetypes/FileHeaderOnly.hxx"
#include "../../lattice/filetypes/FilePlain.hxx"
#include "../../lattice/filetypes/FileVogt.hxx"
#include "../../lattice/filetypes/filetype_typedefs.h"
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;

// lattice setup
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;


typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;


void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);
	s.calculateNeighbourTable( nnt );
}

//
//__global__ void projectSU3( Real* U )
//{
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteCoord<4,FULL_SPLIT> s(size);
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	s.setLatticeIndex( site );
//
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		TLink linkUp( U, s, mu );
//		SU3<TLink> globUp( linkUp );
//
//		globUp.projectSU3();
//	}
//}
//
//
//__global__ void __launch_bounds__(256,4) orStep( Real* U, lat_index_t* nn, bool parity, float orParameter, int counter=0  )
//{
//	typedef GpuLandauPattern< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
//	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;
//
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteIndex<4,true> s(size);
//	s.nn = nn;
//
//	const bool updown = threadIdx.x / 128;
//	const short mu = (threadIdx.x % 128) / 32;
//	const short id = (threadIdx.x % 128) % 32;
//
//	int site = blockIdx.x * blockDim.x/8 + id;
//	if( parity == 1 ) site += s.getLatticeSize()/2;
//
//	s.setLatticeIndex( site );
//	if( updown==1 )
//	{
//		s.setNeighbour(mu,false);
//	}
//
////	if(id == 0) printf("bin in or\n");
//
//	Matrix<complex,Nc> locMat;
//	SU3<Matrix<complex,Nc> > locU(locMat);
//
//	TLinkIndex link( U, s, mu );
//
//	SU3<TLinkIndex> globU( link );
//
//	// make link local
//	locU.assignWithoutThirdLine(globU);
//	locU.reconstructThirdLine();
//
//	// define the update algorithm
//	SrUpdate overrelax( orParameter );
//	GaugeFixingSubgroupStep<SU3<Matrix<complex,Nc> >, SrUpdate, MAG> subgroupStep( &locU, overrelax, id, mu, updown, counter );
//
//	// do the subgroup iteration
//	SU3<Matrix<complex,Nc> >::perSubgroup( subgroupStep );
//
//	// copy link back
//	globU=locU;
//
//	// project back
////	globU.projectSU3withoutThirdRow();
//}




//Real calculatePolyakovLoopAverage( Real *U )
//{
//	Matrix<complex,3> tempMat;
//	SU3<Matrix<complex,3> > temp( tempMat );
//	Matrix<complex,3> temp2Mat;
//	SU3<Matrix<complex,3> > temp2( temp2Mat );
//
//	SiteCoord<Ndim,true> s( HOST_CONSTANTS::SIZE );
//
//	complex result(0,0);
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

	MAGKernelsSU3::initCacheConfig();

	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, options.getDeviceNumber() );
	cudaSetDevice(options.getDeviceNumber());

	printf("\nDevice %d: \"%s\"\n", options.getDeviceNumber(), deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);


	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly( options.getReinterpret() );
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt( options.getReinterpret() );
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain( options.getReinterpret() );


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for configuration
	Real* dU;
	cudaMalloc( &dU, arraySize*sizeof(Real) );

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


	
	GaugeFixingStats<Ndim,Nc,MAGKernelsSU3,AVERAGE> gaugeStats( dU, HOST_CONSTANTS::SIZE );


	// timer to measure kernel times
	Chronotimer kernelTimer;
	kernelTimer.reset();
	kernelTimer.start();

	double saTotalKernelTime = 0;
	double srTotalKernelTime = 0;
	long srTotalStepnumber = 0;
	double orTotalKernelTime = 0; // sum up total kernel time for OR
	long orTotalStepnumber = 0;

	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
	{
		bool loadOk;

		if( !options.isSetHot() ) // load a file
		{
			switch( options.getFType() )
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
		}
		else // or initialize with a hot configuration (ignore file options)
		{
			CommonKernelsSU3::setHot( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize(), options.getSeed(), PhiloxWrapper::getNextCounter() );
		}

		double bestGff = 0.0;
		for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
		{
			// we copy from host in every gaugecopy step to have a cleaner configuration (concerning numerical errors)
			if( !options.isSetHot() ) cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

			if( options.isRandomTrafo() ) // I'm an optimist! This should be called isRandomTrafo()!
			{
				MAGKernelsSU3::randomTrafo(numBlocks,threadsPerBlock,dU, dNn, 0, options.getSeed(), PhiloxWrapper::getNextCounter() );
				MAGKernelsSU3::randomTrafo(numBlocks,threadsPerBlock,dU, dNn, 1, options.getSeed(), PhiloxWrapper::getNextCounter() );
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


				MAGKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 0, temperature, options.getSeed(), PhiloxWrapper::getNextCounter() );
				MAGKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 1, temperature, options.getSeed(), PhiloxWrapper::getNextCounter() );

				for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
				{
					MAGKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 0 );
					MAGKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 1 );
				}


				if( i % options.getCheckPrecision() == 0 )
				{
					CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );

					gaugeStats.generateGaugeQuality();
					printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
				}
				temperature -= tempStep;
			}
			cudaDeviceSynchronize();
			kernelTimer.stop();
			saTotalKernelTime += kernelTimer.getTime();

			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getSrMaxIter(); i++ )
			{

				MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getSrParameter(), options.getSeed(), PhiloxWrapper::getNextCounter() );
				MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getSrParameter(), options.getSeed(), PhiloxWrapper::getNextCounter() );

				if( i % options.getCheckPrecision() == 0 )
				{
					CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
				}
				srTotalStepnumber++;
			}
			cudaDeviceSynchronize();
			kernelTimer.stop();
			srTotalKernelTime += kernelTimer.getTime();

			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{

				MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getOrParameter() );
				MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getOrParameter() );

				if( i % options.getCheckPrecision() == 0 )
				{
					CommonKernelsSU3::projectSU3( numBlocks*2,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
				}

				orTotalStepnumber++;
			}
			cudaDeviceSynchronize();
			kernelTimer.stop();
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
		if( !options.isSetHot() )
		{
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
	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;

	long hbFlops = 2252+86-8;
	long microFlops = 2252+14-8;
	cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

	long srFlops = 2252+32-8;
	cout << "Stochastic Relaxation: " << (double)((long)srFlops*(long)s.getLatticeSize()*(long)srTotalStepnumber)/srTotalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(srTotalStepnumber)*(long)sizeof(Real))/srTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


	long orFlops = 2252+22-8;
	cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

}
