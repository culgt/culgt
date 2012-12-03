/*
 * MultiGPU_MPI_LandauGaugeFixingSU3_4D: main.cpp
 *
 *  Created on: Nov. 16, 2012
 *      Author: vogt&schroeck
 */

#include <iostream>
#include <math.h>
#include <sstream>
#include <cuda_runtime.h>
#include <mpi.h>
#ifndef OSX
#include "malloc.h"
#endif
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/access_pattern/GpuCoulombPatternParity.hxx"
#include "../../lattice/access_pattern/GpuLandauPatternParity.hxx"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/SiteIndex.hxx"
#include "../../lattice/Link.hxx"
#include "../../lattice/SU3.hxx"
#include "../../lattice/Matrix.hxx"
#include "../../lattice/LinkFile.hxx"
#include "../../util/timer/Chronotimer.h"
#include "../../lattice/filetypes/FileHeaderOnly.hxx"
#include "../../lattice/filetypes/FilePlain.hxx"
#include "../../lattice/filetypes/FileVogt.hxx"
#include "../../lattice/filetypes/filetype_typedefs.h"
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include "../../lattice/gaugefixing/GlobalConstants.hxx"
// #include "../../lattice/gaugefixing/LandauKernelsSU3.hxx"
#include "../program_options/ProgramOptions.hxx"
#include "../program_options/FileIterator.hxx"

#include "MultiGPU_MPI_LandauGaugeFixingSU3_4D.h"

// #include "./MPI_ProcInfo.h"
#include "./MultiGPU_MPI_Communicator.hxx"


// // MPI error handling macro
// #define MPI_CHECK( call) \
//     if((call) != MPI_SUCCESS) { \
//         cerr << "MPI error calling \""#call"\"\n"; \
//         MPI_Abort(MPI_COMM_WORLD, (-1) ); }

using namespace std;

//TODO where to put these constants?
// const lat_dim_t Ndim = 4;
// const short Nc = 3;
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
// const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

// lattice setup
const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
const lat_coord_t sizeTimeslice[Ndim] = {1,Nx,Ny,Nz};



typedef GpuCoulombPatternParity<SiteCoord<Ndim,TIMESLICE_SPLIT>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;

void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(size);
	s.calculateNeighbourTable( nnt );
}





int main(int argc, char* argv[])
{
	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;
	
// 	// initialize MPI communication
// 	int nprocs, rank, namelen;
//   char processor_name[MPI_MAX_PROCESSOR_NAME];
// 
//   MPI_CHECK( MPI_Init(&argc, &argv) );
//   MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) );
//   MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
//   MPI_CHECK( MPI_Get_processor_name(processor_name, &namelen) );
// 	printf("Process %d on %s out of %d alive.\n", rank, processor_name, nprocs);
	
	
	//TODO remove later
// 	MPI_Request request1, request2;
// 	MPI_Status  status;

  
	
	// instantiate object of MPI communicator
	MultiGPU_MPI_Communicator comm(argc,argv);
	
	// init. the devices
// 	comm.initDevice( rank%4 );
	
// 		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	
// 	MPI_Finalize();
// 	return 0;
	
	
	// array to distribute the timeslices to the devices
	// TODO put somewhere
// 	int theProcess[Nt];
// 	for( int k=0; k<nprocs; k++ )
// 		for( int t = k*Nt/nprocs; t < (k+1)*Nt/nprocs; t++ )
// 			theProcess[t] = k;
		
		
		
// 	int lRank = ( rank - 1 + nprocs ) % nprocs;
// 	int rRank = ( rank + 1 ) % nprocs;
// 	bool isMaster = ( rank == 0 ? true : false );
// 
// 		
// // 	for( int t=0; t<Nt; t++ )
// // 		if( isMaster) cout << t << ' ' << theProcess[t] << endl;
// 		
// 	
// 	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	
	
// 	// distribute the timeslices over the threads: get tmin and tmax
// 	int tmin = rank*Nt/nprocs;
// 	int tmax = (rank+1)*Nt/nprocs;
// 	int numbSlices = tmax-tmin;
// 	
// 	// 6 intervals for async. mem. copies
// 	int tPart_beg[6];
// 	int tPart_end[6];
// 	
// 	// init.
// 	for( int l=0; l<6; l++ )
// 	{
// 		tPart_beg[l] = tmin+1;
// 		tPart_end[l] = tmin+1;
// 	}
// 	
// 
// 	
// 	for( int t=1; t<numbSlices; t++ )
// 	{
// 		for( int l=0; l<6; l++ )
// 		{
// 			if( l == (t-1)%6 )
// 			{
// 				tPart_end[l] += 1;
// 			}
// 			if( l > (t-1)%6 )
// 			{
// 				tPart_beg[l] += 1;
// 				tPart_end[l] += 1;
// 			}
// 		}
// 	}
// 		
// 	
// 	
// 	printf("Process %d: numbSlices %d\n", rank, numbSlices );
// 	
// 	for( int l=0; l<6; l++ )
// 	{
// 		printf("Process %d: tPart_beg[%d] = %d, tPart_end[%d] = %d\n", rank, l, tPart_beg[l], l, tPart_end[l] );
// 	}
	
//  	printf("Process %d: tmin = %d, tm1 = %d, tm2 = %d, tm3 = %d, tm4 = %d, tm5 = %d, tmax = %d\n", rank, tmin, tm1, tm2, tm3,tm4,  tm5, tmax);
	
// 	MPI_Finalize();
// 	return 0;
	
			
	// how many timeslices does each device have?
// 	int numbTimeSlices=0;
// 	for( int t=tmin; t<tmax; t++ )
// 		numbTimeSlices++;
// 	
// 	printf("Process %d: numbTimeSlices = %d\n", rank, numbTimeSlices);

// 	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
    
	
	// cudaStreams //TODO remove later
// 	cudaStream_t streamStd;
// 	cudaStream_t streamCpy;
// 	cudaStreamCreate( &streamStd );
// 	cudaStreamCreate( &streamCpy );

	Chronotimer allTimer;
	if( comm.isMaster() ) allTimer.reset();

	SiteCoord<4,TIMESLICE_SPLIT> s(HOST_CONSTANTS::SIZE);
	
	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfHeaderOnly;
	lfHeaderOnly.reinterpret = options.getReinterpret();
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfVogt;
	lfVogt.reinterpret = options.getReinterpret();
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfPlain;
	lfPlain.reinterpret = options.getReinterpret();
	
	// allocate Memory
	// host memory for configuration
	Real* U;
	if( comm.isMaster() ) U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for all timeslices
	Real* dU[Nt];
	for( int t=comm.getMinTimeslice(); t<comm.getMaxTimeslice(); t++ )
	{
		cudaMalloc( &dU[t], timesliceArraySize*sizeof(Real) );
	}
	
// 	// halo exchange size
// 	size_t haloSize = timesliceArraySize*sizeof(Real);
// 		
// 	// page-locked host memory for halo timeslices (two per thread)
// 	Real* haloOut[nprocs];
// 	Real* haloIn[nprocs];
//  	cudaHostAlloc( &haloIn[rank],  haloSize, 0 );
// 	cudaHostAlloc( &haloOut[rank], haloSize, 0 );
// 	
// 	// device memory for halo timeslice (one per device)
// 	Real* dHalo[nprocs];
// 	cudaMalloc( &dHalo[rank], haloSize );
	
// 	// device memory for collecting the parts of the gauge fixing functional and divA
// 	double *dGff[nprocs];
// 	double *dA[nprocs];
// 	cudaMalloc( &dGff[rank], s.getLatticeSizeTimeslice()*sizeof(double) );
// 	cudaMalloc( &dA[rank], s.getLatticeSizeTimeslice()*sizeof(double) );

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt[comm.getNumbProcs()];
	cudaMalloc( &dNnt[comm.getRank()], s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	
	// copy neighbour table to device
	cudaMemcpy( dNnt[comm.getRank()], nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t), cudaMemcpyHostToDevice );


	if( comm.isMaster() ) allTimer.start();

// 	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );
	

	

	float totalKernelTime = 0;
	long totalStepNumber = 0;

	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
	{
		// load file
		bool loadOk;
		if( comm.isMaster() )
		{
			cout << "loading " << fi.getFilename() << " as " << options.getFType() << endl;
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
		}
		
		// send (parts of) config. to all processes and transfer to devices
		comm.scatterGaugeField( dU, U );
// 		for( int t=0; t<Nt; t++ )
// 		{
// 			if( comm.isMaster() && theProcess[t] == 0 )
// 			{
// 				cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
// 			}
// 			else if( theProcess[t] == rank )
// 			{
// 				MPI_CHECK( MPI_Recv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status ) );
// 				cudaMemcpy( dU[t], haloIn[rank], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
// 			}
// 			else if( comm.isMaster() )
// 			{
// 				MPI_CHECK( MPI_Send( &U[t*timesliceArraySize], timesliceArraySize, MPI_FLOAT, theProcess[t], 0, MPI_COMM_WORLD ) );
// 			}
// 			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
// 		}
		
// 		// don't load file, fill with random numbers - then projectSU3
// 		if( comm.isMaster() ) cout << "\nWe fill dU with random SU(3) matrices.\n" << endl;
// 		int unsigned counter = 0;//time(NULL);
// 		
// 		for( int t=comm.getMinTimeslice(); t<comm.getMaxTimeslice(); t++ )
// 		{
//  			_set_hot( dU[t], counter ); 
// 			counter++;
// 		}
// 		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		
		
		// calculate and print the gauge quality
		if( comm.isMaster() ) printf( "i:\t\tgff:\t\tdA:\n");
		comm.generateGaugeQuality( dU, dNnt );
		
		//TODO divide somewhere else
		if( comm.isMaster() ) printf( "-\t\t%1.10f\t\t%e\n", comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
// 		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		
			
			
		Chronotimer kernelTimer;
		if( comm.isMaster() ) kernelTimer.reset();
		if( comm.isMaster() ) kernelTimer.start();
		
		// iterate the algorithm
		for( int i = 0; i < options.getOrMaxIter(); i++ )
		{
			for( int evenodd=0; evenodd<2; evenodd++ )
			{
				comm.apply( dU, dNnt, evenodd, OR );
			} 




			// calculate and print the gauge quality
			if( i % options.getCheckPrecision() == 0 )
			{
				comm.generateGaugeQuality( dU, dNnt );
				
				if( comm.isMaster() ) printf( "%d\t\t%1.10f\t\t%e\n", i, comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
// 				MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
			}

			totalStepNumber++;
		}
 		cudaDeviceSynchronize();
		if( comm.isMaster() ) kernelTimer.stop();
		if( comm.isMaster() ) cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		

		
		// send back all timeslices to master
		comm.collectGaugeField( dU, U );
// 		for( int t=0; t<Nt; t++ )
// 		{
// 			if( comm.isMaster() && theProcess[t] == 0 )
// 			{
// 				cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
// 			}
// 			else if( theProcess[t] == rank )
// 			{
// 				MPI_CHECK( MPI_Send( haloIn[rank],  timesliceArraySize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD ) );
// 				cudaMemcpy( haloIn[rank], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
// 			}
// 			else if( comm.isMaster() )
// 			{
// 				MPI_CHECK( MPI_Recv( &U[t*timesliceArraySize], timesliceArraySize, MPI_FLOAT, theProcess[t], 0, MPI_COMM_WORLD, &status ) );
// 			}
// 			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
// 		}
// 		
// 		//saving file
// 		if( comm.isMaster() )
// 		{
// 			cout << "saving " << fi.getOutputFilename() << " as " << options.getFType() << endl;
// 			switch( options.getFType() )
// 			{
// 			case VOGT:
// 				loadOk = lfVogt.save( s, fi.getOutputFilename(), U );
// 				break;
// 			case PLAIN:
// 				loadOk = lfPlain.save( s, fi.getOutputFilename(), U );
// 				break;
// 			case HEADERONLY:
// 				loadOk = lfHeaderOnly.save( s, fi.getOutputFilename(), U );
// 				break;
// 			default:
// 				cout << "Filetype not set to a known value. Exiting";
// 				exit(1);
// 			}
// 		}
	} // end fileIterator

	if( comm.isMaster() ) allTimer.stop();
	if( comm.isMaster() ) cout << "total time: " << allTimer.getTime() << " s" << endl;
	if( comm.isMaster() ) cout << "total kernel time: " << totalKernelTime << " s" << endl;
	if( comm.isMaster() ) cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


// 	MPI_CHECK( MPI_Finalize() );
	return 0;

}
