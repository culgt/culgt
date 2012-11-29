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
#include "../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
#include "../lattice/gaugefixing/GaugeFixingStatsV2.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../lattice/gaugefixing/overrelaxation/SrUpdate.hxx"
#include "../lattice/gaugefixing/MAGKernelsSU3.hxx"
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
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include "../lattice/gaugefixing/GlobalConstants.hxx"
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;



// boost program options setup
//boost::program_options::variables_map options_vm;
//boost::program_options::options_description options_desc("Allowed options");
//
//// parameters from command line or config file
//int nconf;
//int devicenumb;
//long seed; // TODO check datatype
//int orMaxIter;
//int orCheckPrec;
//float orParameter;
//float orPrecision;
//int saSteps;
//float saMin;
//float saMax;
//int gaugeCopies;
//string fileEnding;
//string postFixLabel;
//string fileBasename;
//int fileStartnumber;
//int fileStepsize;
//int fileNumberformat;
//string configFile;
//bool noRandomTrafo;
//FileType fileType;


// lattice setup
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;


typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;


//__device__ inline Real cuFabs( Real a )
//{
//	return (a>0)?(a):(-a);
//}

void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);
	s.calculateNeighbourTable( nnt );
}


__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteCoord<4,FULL_SPLIT> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );

		globUp.projectSU3();
	}
}
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

	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


	// TODO maybe we should choose the filetype on compile time
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

	
	// instantiate GaugeFixingStats object
//	lat_coord_t *devicePointerToSize;
//	cudaGetSymbolAddress( (void**)&devicePointerToSize, "dSize" );
	GaugeFixingStats<Ndim,Nc,MAGKernelsSU3,AVERAGE> gaugeStats( dU, HOST_CONSTANTS::SIZE );


	double totalKernelTime = 0;

	long totalStepNumber = 0;

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

		// copying configuration ...
		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

		// calculate and print the gauge quality
		printf( "i:\t\tgff:\t\tdA:\n");
		gaugeStats.generateGaugeQuality();
		printf( "   \t\t%1.10f\t\t%e\n", gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

		Chronotimer kernelTimer;
		kernelTimer.reset();
		kernelTimer.start();



		float temperature = options.getSaMax();
		float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

		for( int i = 0; i < options.getSaSteps(); i++ )
		{


			MAGKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 0, temperature, PhiloxWrapper::getNextCounter() );
			MAGKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 1, temperature, PhiloxWrapper::getNextCounter() );

			for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
			{
				MAGKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 0 );
				MAGKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 1 );
			}


			if( i % options.getCheckPrecision() == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );


				gaugeStats.generateGaugeQuality();
//				CudaError::getLastError( "generateGaugeQuality error" );
				printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
			}
			temperature -= tempStep;
		}

		for( int i = 0; i < options.getSrMaxIter(); i++ )
		{

			MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getSrParameter(), PhiloxWrapper::getNextCounter() );
			MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getSrParameter(), PhiloxWrapper::getNextCounter() );

			if( i % options.getCheckPrecision() == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );
				gaugeStats.generateGaugeQuality();
				printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

				if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
			}
		}

		for( int i = 0; i < options.getOrMaxIter(); i++ )
		{

			MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getOrParameter() );
			MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getOrParameter() );

			if( i % options.getCheckPrecision() == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );
				gaugeStats.generateGaugeQuality();
				printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

				if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
			}

			totalStepNumber++;
		}

		cudaThreadSynchronize();
		kernelTimer.stop();
		cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();

		// copy back
		cudaMemcpy( U, dU, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		
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
	cout << "total kernel time: " << totalKernelTime << " s" << endl;

	cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

}
