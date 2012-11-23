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


// MPI error handling macro
#define MPI_CHECK( call) \
    if((call) != MPI_SUCCESS) { \
        cerr << "MPI error calling \""#call"\"\n"; \
        MPI_Abort(MPI_COMM_WORLD, (-1) ); }

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;


// lattice setup
const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
const lat_coord_t sizeTimeslice[Ndim] = {1,Nx,Ny,Nz};

const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

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
	
	// initialize MPI communication
	int nprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_CHECK( MPI_Init(&argc, &argv) );
  MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) );
  MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
  MPI_CHECK( MPI_Get_processor_name(processor_name, &namelen) );
	
	MPI_Request request1, request2;
	MPI_Status  status;

  printf("Process %d on %s out of %d alive.\n", rank, processor_name, nprocs);

	int lRank = ( rank - 1 + nprocs ) % nprocs;
	int rRank = ( rank + 1 ) % nprocs;
	bool isMaster = ( rank == 0 ? true : false );
	
	// array to distribute the timeslices to the devices
	int theProcess[Nt];
	for( int k=0; k<nprocs; k++ )
		for( int t = k*Nt/nprocs; t < (k+1)*Nt/nprocs; t++ )
			theProcess[t] = k;
		
// 	for( int t=0; t<Nt; t++ )
// 		if( isMaster) cout << t << ' ' << theProcess[t] << endl;
		

	initDevice( rank );
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	
	
	// distribute the timeslices over the threads: get tmin and tmax
	int tmin = rank*Nt/nprocs;
	int tmax = (rank+1)*Nt/nprocs;
	int numbSlices = tmax-tmin;
	
	// 6 intervals for async. mem. copies
	int tPart_beg[6];
	int tPart_end[6];
	
	// init.
	for( int l=0; l<6; l++ )
	{
		tPart_beg[l] = tmin+1;
		tPart_end[l] = tmin+1;
	}
	

	
	for( int t=1; t<numbSlices; t++ )
	{
		for( int l=0; l<6; l++ )
		{
			if( l == (t-1)%6 )
			{
				tPart_end[l] += 1;
			}
			if( l > (t-1)%6 )
			{
				tPart_beg[l] += 1;
				tPart_end[l] += 1;
			}
		}
	}
		
	
	
	printf("Process %d: numbSlices %d\n", rank, numbSlices );
	
	for( int l=0; l<6; l++ )
	{
		printf("Process %d: tPart_beg[%d] = %d, tPart_end[%d] = %d\n", rank, l, tPart_beg[l], l, tPart_end[l] );
	}
	
//  	printf("Process %d: tmin = %d, tm1 = %d, tm2 = %d, tm3 = %d, tm4 = %d, tm5 = %d, tmax = %d\n", rank, tmin, tm1, tm2, tm3,tm4,  tm5, tmax);
	
// 	MPI_Finalize();
// 	return 0;
	
			
	// how many timeslices does each device have?
	int numbTimeSlices=0;
	for( int t=tmin; t<tmax; t++ )
		numbTimeSlices++;
	
	printf("Process %d: numbTimeSlices = %d\n", rank, numbTimeSlices);
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
    
	
	// cudaStreams
	cudaStream_t streamStd;
	cudaStream_t streamCpy;
	cudaStreamCreate( &streamStd );
	cudaStreamCreate( &streamCpy );

	Chronotimer allTimer;
	if( isMaster ) allTimer.reset();

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
	if( isMaster ) U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for all timeslices
	Real* dU[Nt];
	for( int t=tmin; t<tmax; t++ )
	{
		cudaMalloc( &dU[t], timesliceArraySize*sizeof(Real) );
	}
	
	// halo exchange size
	size_t haloSize = timesliceArraySize*sizeof(Real);
		
	// page-locked host memory for halo timeslices (two per thread)
	Real* haloOut[nprocs];
	Real* haloIn[nprocs];
 	cudaHostAlloc( &haloIn[rank],  haloSize, 0 );
	cudaHostAlloc( &haloOut[rank], haloSize, 0 );
	
	// device memory for halo timeslice (one per device)
	Real* dHalo[nprocs];
	cudaMalloc( &dHalo[rank], haloSize );
	
	// device memory for collecting the parts of the gauge fixing functional and divA
	double *dGff[nprocs];
	double *dA[nprocs];
	cudaMalloc( &dGff[rank], s.getLatticeSizeTimeslice()*sizeof(double) );
	cudaMalloc( &dA[rank], s.getLatticeSizeTimeslice()*sizeof(double) );

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt[nprocs];
	cudaMalloc( &dNnt[rank], s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	
	// copy neighbour table to device
	cudaMemcpy( dNnt[rank], nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t), cudaMemcpyHostToDevice );


	if( isMaster ) allTimer.start();

// 	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );

	float totalKernelTime = 0;
	long totalStepNumber = 0;

	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
	{
		// load file
		if( isMaster )
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
		}
		
		// send (parts of) config. to all processes and transfer to devices
		for( int t=0; t<Nt; t++ )
		{
			if( isMaster && theProcess[t] == 0 )
			{
				cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
			}
			else if( theProcess[t] == rank )
			{
				MPI_CHECK( MPI_Recv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status ) );
				cudaMemcpy( dU[t], haloIn[rank], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
			}
			else if( isMaster )
			{
				MPI_CHECK( MPI_Send( &U[t*timesliceArraySize], timesliceArraySize, MPI_FLOAT, theProcess[t], 0, MPI_COMM_WORLD ) );
			}
			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		}
		
// 		// don't load file, fill with random numbers - then projectSU3
// 		if( isMaster ) cout << "\nWe fill dU with random SU(3) matrices.\n" << endl;
// 		int unsigned counter = 0;//time(NULL);
// 		
// 		for( int t=tmin; t<tmax; t++ )
// 		{
//  			_set_hot( dU[t], counter ); 
// 			counter++;
// 		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		
		
		// calculate and print the gauge quality
		if( isMaster ) printf( "i:\t\tgff:\t\tdA:\n");
		// host memory (two fields to collect the results)
		double gff[2], A[2];
		double gffRes, ARes;
		

		gff[1]=0.0; 
		A[1]  =0.0;
	
		//exchange halos
		cudaMemcpy( haloOut[rank], dU[tmax-1], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		cudaDeviceSynchronize();
		MPI_CHECK( MPI_Irecv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2) );		
		MPI_CHECK( MPI_Isend( haloOut[rank], timesliceArraySize, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1) );
		MPI_CHECK( MPI_Wait( &request1, &status ) );
		MPI_CHECK( MPI_Wait( &request2, &status ) );
		cudaMemcpy( dHalo[rank], haloIn[rank], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		
		// call kernel wrapper for all t
		for( int t=tmin; t<tmax; t++ )
		{
			if( t == tmin )
				_generateGaugeQualityPerSite( dU[tmin], dHalo[rank], dNnt[rank], dGff[rank], dA[rank] );
			else
				_generateGaugeQualityPerSite( dU[t], dU[t-1], dNnt[rank], dGff[rank], dA[rank] );
			
			_averageGaugeQuality( dGff[rank], dA[rank] );
			cudaMemcpy( &gff[0], dGff[rank], sizeof(double), cudaMemcpyDeviceToHost );
			cudaMemcpy( &A[0], dA[rank],     sizeof(double), cudaMemcpyDeviceToHost );
			gff[1]+=gff[0];
			A[1]  +=A[0];
		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		
		MPI_CHECK( MPI_Reduce( &gff[1], &gffRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
		MPI_CHECK( MPI_Reduce( &A[1],   &ARes,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
		
		gffRes /= (double)Nt;
		ARes   /= (double)Nt;
		
		if( isMaster ) printf( "-\t\t%1.10f\t\t%e\n", gffRes, ARes );
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		
			
			
		Chronotimer kernelTimer;
		if( isMaster ) kernelTimer.reset();
		if( isMaster ) kernelTimer.start();
		
		// iterate the algorithm
		for( int i = 0; i < options.getOrMaxIter(); i++ )
		{
			for( int parity=0; parity<2; parity++ )
			{
				if( nprocs > 1 )
				{
					int p_offset = parity?timesliceArraySize/2:0;
					
					// halo exchange forward step 1
					cudaMemcpyAsync( haloOut[rank]+p_offset, dU[tmax-1]+p_offset, haloSize/12, cudaMemcpyDeviceToHost, streamCpy );
					for( int t = tPart_beg[2]; t < tPart_end[2]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
					cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
					
					// halo exchange forward step 2
					MPI_CHECK( MPI_Irecv( haloIn[rank]+p_offset,  timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2) );	
					MPI_CHECK( MPI_Isend( haloOut[rank]+p_offset, timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1) );
					for( int t = tPart_beg[0]; t < tPart_end[0]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
					MPI_CHECK( MPI_Wait( &request1, &status ) );
					MPI_CHECK( MPI_Wait( &request2, &status ) );
					
					// halo exchange forward step 3
					cudaMemcpyAsync( dHalo[rank]+p_offset, haloIn[rank]+p_offset, haloSize/12, cudaMemcpyHostToDevice, streamCpy );
					for( int t = tPart_beg[3]; t < tPart_end[3]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
					
					// now call kernel wrapper for tmin with dU[t-1] replaced by dHalo[rank]
					_orStep( dU[tmin], dHalo[rank], dNnt[rank], parity, options.getOrParameter(), streamCpy );
					
					// halo exchange back step 1
					cudaMemcpyAsync( haloOut[rank]+p_offset, dHalo[rank]+p_offset, haloSize/12, cudaMemcpyDeviceToHost, streamCpy );
					for( int t = tPart_beg[4]; t < tPart_end[4]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
					cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
					
					// halo exchange back step 2
					MPI_CHECK( MPI_Irecv( haloIn[rank]+p_offset,  timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request2) );	
					MPI_CHECK( MPI_Isend( haloOut[rank]+p_offset, timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request1) );
					for( int t = tPart_beg[1]; t < tPart_end[1]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
					MPI_CHECK( MPI_Wait( &request1, &status ) );
					MPI_CHECK( MPI_Wait( &request2, &status ) );
					
					// halo exchange back step 3
					cudaMemcpyAsync( dU[tmax-1]+p_offset, haloIn[rank]+p_offset, haloSize/12, cudaMemcpyHostToDevice, streamCpy );
					for( int t = tPart_beg[5]; t < tPart_end[5]; t++ )
					{
						_orStep( dU[t], dU[t-1], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}		
					cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
					
					MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
				}
				else // nproc == 1
				{
					for( int t=tmin; t<tmax; t++ )
					{
						int tDw = ( t > 0 )?( t - 1 ):( Nt - 1 );
						_orStep( dU[t], dU[tDw], dNnt[rank], parity, options.getOrParameter(), streamStd );
					}
				} // end if nproc > 1
			} // end for parity




			// calculate and print the gauge quality
			if( i % options.getOrCheckPrecision() == 0 )
			{
				gff[1]=0.0; 
				A[1]  =0.0;
			
				//exchange halos
				cudaMemcpy( haloOut[rank], dU[tmax-1], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
				cudaDeviceSynchronize();
				MPI_CHECK( MPI_Irecv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2) );		
				MPI_CHECK( MPI_Isend( haloOut[rank], timesliceArraySize, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1) );
				MPI_CHECK( MPI_Wait( &request1, &status ) );
				MPI_CHECK( MPI_Wait( &request2, &status ) );
				cudaMemcpy( dHalo[rank], haloIn[rank], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
				
				// call kernel wrapper for all t
				for( int t=tmin; t<tmax; t++ )
				{
					if( t == tmin )
						_generateGaugeQualityPerSite( dU[tmin], dHalo[rank], dNnt[rank], dGff[rank], dA[rank] );
					else
						_generateGaugeQualityPerSite( dU[t], dU[t-1], dNnt[rank], dGff[rank], dA[rank] );
					
					_averageGaugeQuality( dGff[rank], dA[rank] );
					cudaMemcpy( &gff[0], dGff[rank], sizeof(double), cudaMemcpyDeviceToHost );
					cudaMemcpy( &A[0], dA[rank],     sizeof(double), cudaMemcpyDeviceToHost );
					gff[1]+=gff[0];
					A[1]  +=A[0];
				}
				MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
				
				MPI_CHECK( MPI_Reduce( &gff[1], &gffRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
				MPI_CHECK( MPI_Reduce( &A[1],   &ARes,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
				
				gffRes /= (double)Nt;
				ARes   /= (double)Nt;
				
				if( isMaster ) printf( "%d\t\t%1.10f\t\t%e\n", i, gffRes, ARes );
				MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
			}

			totalStepNumber++;
		}
 		cudaDeviceSynchronize();
		if( isMaster ) kernelTimer.stop();
		if( isMaster ) cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		

		// copy back all timeslices
// 		for( int t=0; t<Nt; t++ )
// 			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
// 		// copy back
// 		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );
// 		
// 		//saving file
// 		cout << "saving " << fi.getOutputFilename() << " as " << options.getFType() << endl;
// 
// 		switch( options.getFType() )
// 		{
// 		case VOGT:
// 			loadOk = lfVogt.save( s, fi.getOutputFilename(), U );
// 			break;
// 		case PLAIN:
// 			loadOk = lfPlain.save( s, fi.getOutputFilename(), U );
// 			break;
// 		case HEADERONLY:
// 			loadOk = lfHeaderOnly.save( s, fi.getOutputFilename(), U );
// 			break;
// 		default:
// 			cout << "Filetype not set to a known value. Exiting";
// 			exit(1);
// 		}
	}

	if( isMaster ) allTimer.stop();
	if( isMaster ) cout << "total time: " << allTimer.getTime() << " s" << endl;
	if( isMaster ) cout << "total kernel time: " << totalKernelTime << " s" << endl;
	if( isMaster ) cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


	MPI_CHECK( MPI_Finalize() );
	return 0;

}
