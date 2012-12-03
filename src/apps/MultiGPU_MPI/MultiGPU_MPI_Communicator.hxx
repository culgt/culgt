/*
 *  MultiGPU_MPI_Communicator.hxx
 *
 *  Created on: Nov 29, 2012
 *      Author: schroeck
 * 
 *  compile with mpicc
 * 
 */

#ifndef MULTIGPU_MPI_COMMUNICATOR_HXX_
#define MULTIGPU_MPI_COMMUNICATOR_HXX_

#include "../../util/datatype/datatypes.h"
#include "../../util/datatype/lattice_typedefs.h"
#include "../../lattice/gaugefixing/GlobalConstants.hxx"
#include "./MPI_ProcInfo.h"
#include "./MultiGPU_MPI_LandauKernelsSU3.h"
#include "./MultiGPU_MPI_Reduce.h"
#include <stdio.h>

// MPI error handling macro
#define MPI_CHECK( call) \
    if((call) != MPI_SUCCESS) { \
        cerr << "MPI error calling \""#call"\"\n"; \
        MPI_Abort(MPI_COMM_WORLD, (-1) ); }
        

//TODO where to put these constants?
const int Ndim = 4;
const int Nc = 3;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
const size_t timesliceSize = timesliceArraySize*sizeof(Real); 
        

class MultiGPU_MPI_Communicator
{
public:
	// constructor
	MultiGPU_MPI_Communicator( int argc, char** argv );
	// destructor
	~MultiGPU_MPI_Communicator();
	// scatter the gauge field from 'master' to all other processes
	void scatterGaugeField( Real **dU, Real *U );
	// collect the gauge field from all processes to 'master'
	void collectGaugeField( Real **dU, Real *U );
	// get number of processes in MPI universe
	int getNumbProcs();
	// get rank
	int getRank();
	// get left neighbor rank
	int getLeft();
	// get right neighbor rank
	int getRight();
	// get min. timeslice the process takes care of
	int getMinTimeslice();
	// get max. timeslice the process takes care of
	int getMaxTimeslice();
	// is the process master?
	bool isMaster();
	// get numb. of timeslices the process takes care of
	int getNumbTimeslices();
	// get first slice of the six parts to hide the 6 parts of comm.
	int getStartPart( int beg );
	// get last slice of the six parts to hide the 6 parts of comm.
	int getEndPart( int end );
	// apply: applies algorithm, takes care of communication between devices
	void apply( Real** dU, lat_index_t** dNnt, bool evenodd, enum AlgoType algorithm );
	// generate the gauge quality
	void generateGaugeQuality( Real** dU, lat_index_t** dNnt );
	// get the current value of the gauge functional
	double getCurrentGff();
	// get the current value of the gauge precision
	double getCurrentA();
	
private:
	// assign a device to the process
	void initDevice( const int device );
	int nprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	bool master;
	int lRank;
	int rRank;
	
	// MPI comm.
	MPI_Request request1, request2;
	MPI_Status  status;
	
	// useful variables
	int tmin;
	int tmax;
	int numbSlices;
	int startPart[6];
	int endPart[6];
	int theProcess[Nt];
	
	// device memory for collecting the parts of the gauge fixing functional and divA
	double *dGff;
	double *dA;
	double currentGff;
	double currentA;
	
	// halos
	Real* haloOut;
	Real* haloIn;
	Real* dHalo;
	
	// cudaStreams 
	cudaStream_t streamStd;
	cudaStream_t streamCpy;
};

MultiGPU_MPI_Communicator::MultiGPU_MPI_Communicator( int argc, char** argv )
{
	
	// initialize MPI communication
  MPI_CHECK( MPI_Init(&argc, &argv) );
  MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) );
  MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
  MPI_CHECK( MPI_Get_processor_name(processor_name, &namelen) );
	printf("Process %d on %s out of %d alive.\n", rank, processor_name, nprocs);
	
	// init. some variables
	lRank = ( rank - 1 + nprocs ) % nprocs;
	rRank = ( rank + 1 ) % nprocs;
	master = ( rank == 0 ? true : false );
	tmin = rank*Nt/nprocs;
	tmax = (rank+1)*Nt/nprocs;
	numbSlices = tmax-tmin;
	
	// set boarders of six parts to hide the comm. (apply,generateGaugeQuality)
	for( int l=0; l<6; l++ )
	{
		startPart[l] = tmin+1;
		endPart[l] = tmin+1;
	}
	for( int t=1; t<numbSlices; t++ )
	{
		for( int l=0; l<6; l++ )
		{
			if( l == (t-1)%6 )
			{
				endPart[l] += 1;
			}
			if( l > (t-1)%6 )
			{
				startPart[l] += 1;
				endPart[l] += 1;
			}
		}
	}
	
	// to which process belongs timeslice 't':
	for( int k=0; k<nprocs; k++ )
		for( int t = k*Nt/nprocs; t < (k+1)*Nt/nprocs; t++ )
			theProcess[t] = k;
	

	// print some information
	printf("Process %d: numbSlices %d\n", rank, numbSlices );
	for( int l=0; l<6; l++ )
	{
		printf("Process %d: startPart[%d] = %d, endPart[%d] = %d\n", rank, l, startPart[l], l, endPart[l] );
	}
	
	// init. the device
	initDevice( (rank%4)+1 );
	
	// init. cuda streams
	cudaStreamCreate( &streamStd );
	cudaStreamCreate( &streamCpy );
	
	// allocate memory for halo exchange
	//TODO timesliceArraySize etc. somewhere else
// 	static const int Ndim = 4;
// 	static const int Nc = 3;
// 	static const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 	size_t timesliceSize = timesliceArraySize*sizeof(Real);
		
	// page-locked host memory for halo timeslices (two per thread)
 	cudaHostAlloc( &haloIn,  timesliceSize, 0 );
	cudaHostAlloc( &haloOut, timesliceSize, 0 );
	
	// device memory for halo timeslice (one per device)
	cudaMalloc( &dHalo, timesliceSize );
	
	// device memory for gauge quality
	cudaMalloc( &dGff, Nx*Ny*Nz*sizeof(double)/2 );
	cudaMalloc( &dA,   Nx*Ny*Nz*sizeof(double)/2 );
	
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}

MultiGPU_MPI_Communicator::~MultiGPU_MPI_Communicator()
{
	// free halo memory
 	cudaFreeHost( haloIn );
	cudaFreeHost( haloOut );
	cudaFree( dHalo );
	
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	MPI_CHECK( MPI_Finalize() );
}

void MultiGPU_MPI_Communicator::initDevice( const int device )
{
	cudaSetDevice(device);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);
	printf("\nDevice %d: \"%s\"\n", device, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);
}

void MultiGPU_MPI_Communicator::scatterGaugeField( Real **dU, Real *U )
{
// 	static const int Ndim = 4;
// 	static const int Nc = 3;
// 	static const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 	static const size_t timesliceSize = timesliceArraySize*sizeof(Real);
	
	//TODO clear this
	// page-locked host memory for halo timeslices (two per thread)
// 	static Real* haloOut[32];
// 	static Real* haloIn[32];
	
	// device memory for halo timeslice (one per device)
// 	static Real* dHalo[32];
	
	
	static const int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

	// MPI comm.
// 	static MPI_Request request1, request2;
// 	static MPI_Status  status;

	
// 	cudaHostAlloc( &haloIn,  timesliceSize, 0 );
// 	cudaHostAlloc( &haloOut, timesliceSize, 0 );
// 	cudaMalloc( &dHalo, timesliceSize );
	
	
	for( int t=0; t<Nt; t++ )
	{
		if( master && theProcess[t] == 0 )
		{
			cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		}
		else if( theProcess[t] == rank )
		{
			MPI_CHECK( MPI_Recv( haloIn,  timesliceArraySize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status ) );
			cudaMemcpy( dU[t], haloIn, timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		}
		else if( master )
		{
			MPI_CHECK( MPI_Send( &U[t*timesliceArraySize], timesliceArraySize, MPI_FLOAT, theProcess[t], 0, MPI_COMM_WORLD ) );
		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	}
}

void MultiGPU_MPI_Communicator::collectGaugeField( Real **dU, Real *U )
{
// 	static const int Ndim = 4;
// 	static const int Nc = 3;
// 	static const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 	static const size_t timesliceSize = timesliceArraySize*sizeof(Real);
	
	//TODO clear this
	// page-locked host memory for halo timeslices (two per thread)
// 	static Real* haloOut[32];
// 	static Real* haloIn[32];
	
	// device memory for halo timeslice (one per device)
// 	static Real* dHalo[32];
	
	
	static const int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

// 	// MPI comm.
// 	static MPI_Request request1, request2;
// 	static MPI_Status  status;

	
// 	cudaHostAlloc( &haloIn,  timesliceSize, 0 );
// 	cudaHostAlloc( &haloOut, timesliceSize, 0 );
// 	cudaMalloc( &dHalo, timesliceSize );
	
	
	// send back all timeslices to master
	for( int t=0; t<Nt; t++ )
	{
		if( master && theProcess[t] == 0 )
		{
			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		}
		else if( theProcess[t] == rank )
		{
			MPI_CHECK( MPI_Send( haloIn,  timesliceArraySize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD ) );
			cudaMemcpy( haloIn, dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		}
		else if( master )
		{
			MPI_CHECK( MPI_Recv( &U[t*timesliceArraySize], timesliceArraySize, MPI_FLOAT, theProcess[t], 0, MPI_COMM_WORLD, &status ) );
		}
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		}
}

int MultiGPU_MPI_Communicator::getNumbProcs()
{
	return nprocs;
}

int MultiGPU_MPI_Communicator::getRank()
{
	return rank;
}

int MultiGPU_MPI_Communicator::getLeft()
{
	return lRank;
}

int MultiGPU_MPI_Communicator::getRight()
{
	return rRank;
}

int MultiGPU_MPI_Communicator::getMinTimeslice()
{
	return tmin;
}

int MultiGPU_MPI_Communicator::getMaxTimeslice()
{
	return tmax;
}

bool MultiGPU_MPI_Communicator::isMaster()
{
	return master;
}

int MultiGPU_MPI_Communicator::getNumbTimeslices()
{
	return numbSlices;
}

int MultiGPU_MPI_Communicator::getStartPart( int beg )
{
	return startPart[ beg ];
}

int MultiGPU_MPI_Communicator::getEndPart( int end )
{
	return endPart[ end ];
}

inline void MultiGPU_MPI_Communicator::apply( Real** dU, lat_index_t** dNnt, bool evenodd, enum AlgoType algorithm )
{
// 	static const int Ndim = 4;
// 	static const int Nc = 3;
// 	static const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 	static const size_t timesliceSize = timesliceArraySize*sizeof(Real);
	
	static const int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

	// MPI comm.
// 	static MPI_Request request1, request2;
// 	static MPI_Status  status;

// 	// cudaStreams 
// 	static cudaStream_t streamStd;
// 	static cudaStream_t streamCpy;
	
	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_LandauKernelsSU3 oneslicehitter;
	
// 	if( init )
// 	{
// 		cudaStreamCreate( &streamStd );
// 		cudaStreamCreate( &streamCpy );
// // 		oneslicehitter.initCacheConfig();
// 	}
	

	
	if( nprocs > 1 )
	{
		int p_offset = evenodd?timesliceArraySize/2:0;
		
		// halo exchange forward step 1
		cudaMemcpyAsync( haloOut+p_offset, dU[tmax-1]+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
		for( int t = startPart[2]; t < endPart[2]; t++ )
		{
			// call wrapper for one timelice
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );
		}
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		// halo exchange forward step 2
		MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2) );	
		MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1) );
		for( int t = startPart[0]; t < endPart[0]; t++ )
		{
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );
		}
		MPI_CHECK( MPI_Wait( &request1, &status ) );
		MPI_CHECK( MPI_Wait( &request2, &status ) );
		
		// halo exchange forward step 3
		cudaMemcpyAsync( dHalo+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
		for( int t = startPart[3]; t < endPart[3]; t++ )
		{
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );		
		}
		
		// now call kernel wrapper for tmin with dU[t-1] replaced by dHalo
		oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamCpy, dU[tmin], dHalo, dNnt[rank], evenodd ^ (tmin%2) , algorithm );
		
		// halo exchange back step 1
		cudaMemcpyAsync( haloOut+p_offset, dHalo+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
		for( int t = startPart[4]; t < endPart[4]; t++ )
		{
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );
		}
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		// halo exchange back step 2
		MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request2) );	
		MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request1) );
		for( int t = startPart[1]; t < endPart[1]; t++ )
		{
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );
		}
		MPI_CHECK( MPI_Wait( &request1, &status ) );
		MPI_CHECK( MPI_Wait( &request2, &status ) );
		
		// halo exchange back step 3
		cudaMemcpyAsync( dU[tmax-1]+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
		for( int t = startPart[5]; t < endPart[5]; t++ )
		{
			oneslicehitter.applyOneTimeslice( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , algorithm );
		}		
		cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
		
		MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	}
	else // nproc == 1
	{
		for( int t=tmin; t<tmax; t++ )
		{
			int tDw = ( t > 0 )?( t - 1 ):( Nt - 1 );
			oneslicehitter.applyOneTimeslice( numBlocks, threadsPerBlock, streamStd, dU[t], dU[tDw], dNnt[rank], evenodd ^ (t%2), algorithm );
		}
	} // end if nproc > 1
}


inline void MultiGPU_MPI_Communicator::generateGaugeQuality( Real** dU, lat_index_t** dNnt )
{
// 	static const int Ndim = 4;
// 	static const int Nc = 3;
// 	static const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;
// 	static const size_t timesliceSize = timesliceArraySize*sizeof(Real);
	
	static const int threadsPerBlock = NSB; // NSB sites are updated within a block (8 threads are needed per site)
	static const int numBlocks = Nx*Ny*Nz/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

// 	// MPI comm.
// 	static MPI_Request request1, request2;
// 	static MPI_Status  status;
// 
// 	// cudaStreams 
// 	static cudaStream_t streamStd;
// 	static cudaStream_t streamCpy;
	
	// instantiate object of kernel wrapper class
	static MultiGPU_MPI_LandauKernelsSU3 oneslicehitter;
	
	// device memory for collecting the parts of the gauge fixing functional and divA
// 	static double *dGff;
// 	static double *dA;
	
	// reduce to collect dGff, dA
	static Reduce reduce( Nx*Ny*Nz/2 );
	
// 	if( init )
// 	{
// // 		cudaStreamCreate( &streamStd );
// // 		cudaStreamCreate( &streamCpy );
// // 		oneslicehitter.initCacheConfig();
// 		// TODO make half the size (parity)
// 		cudaMalloc( &dGff, Nx*Ny*Nz*sizeof(double)/2 );
// 		cudaMalloc( &dA,   Nx*Ny*Nz*sizeof(double)/2 );
// 	}
	
	double tempGff = 0.0;
	double tempA   = 0.0;
	
	if( nprocs > 1 )
	{
		for( int evenodd=0; evenodd<2; evenodd++ )
		{
			int p_offset = evenodd?timesliceArraySize/2:0;
			
			// halo exchange forward step 1
			cudaMemcpyAsync( haloOut+p_offset, dU[tmax-1]+p_offset, timesliceSize/12, cudaMemcpyDeviceToHost, streamCpy );
			for( int t = startPart[0]; t < endPart[1]; t++ )
			{
				// call wrapper for one timelice
				oneslicehitter.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			cudaDeviceSynchronize(); // to ensure cudaMemcpyAsync finished
			
			// halo exchange forward step 2
			MPI_CHECK( MPI_Irecv( haloIn+p_offset,  timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2) );	
			MPI_CHECK( MPI_Isend( haloOut+p_offset, timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1) );
			for( int t = startPart[2]; t < endPart[3]; t++ )
			{
				oneslicehitter.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			MPI_CHECK( MPI_Wait( &request1, &status ) );
			MPI_CHECK( MPI_Wait( &request2, &status ) );
			
			// halo exchange forward step 3
			cudaMemcpyAsync( dHalo+p_offset, haloIn+p_offset, timesliceSize/12, cudaMemcpyHostToDevice, streamCpy );
			for( int t = startPart[4]; t < endPart[5]; t++ )
			{
				oneslicehitter.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamStd, dU[t], dU[t-1], dNnt[rank], evenodd ^ (t%2) , dGff, dA );		
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
			
			// now call kernel wrapper for tmin with dU[t-1] replaced by dHalo
			oneslicehitter.generateGaugeQualityPerSite( numBlocks,threadsPerBlock, streamCpy, dU[tmin], dHalo, dNnt[rank], evenodd ^ (tmin%2) , dGff, dA );
			tempGff += reduce.getReducedValue( streamCpy, dGff );
			tempA   += reduce.getReducedValue( streamCpy, dA );
				
			cudaDeviceSynchronize();
			MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
		}
	}
	else // nproc == 1
	{
		for( int evenodd=0; evenodd<2; evenodd++ )
			for( int t=tmin; t<tmax; t++ )
			{
				int tDw = ( t > tmin )?( t - 1 ):( tmax - 1 );
				oneslicehitter.generateGaugeQualityPerSite( numBlocks, threadsPerBlock, streamStd, dU[t], dU[tDw], dNnt[rank], evenodd ^ (t%2), dGff, dA );
				tempGff += reduce.getReducedValue( streamStd, dGff );
				tempA   += reduce.getReducedValue( streamStd, dA );
			}
	} // end if nproc > 1
	

	// collect dGff, dA from all devices
	if( nprocs > 1 )
	{
		MPI_CHECK( MPI_Reduce( &tempGff, &currentGff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
		MPI_CHECK( MPI_Reduce( &tempA,   &currentA,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	}
	else
	{
		currentGff = tempGff;
		currentA = tempA;
	}
	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
}

double MultiGPU_MPI_Communicator::getCurrentGff()
{
	return currentGff;
}

double MultiGPU_MPI_Communicator::getCurrentA()
{
	return currentA;
}



#endif /* MULTIGPU_MPI_COMMUNICATOR_HXX_ */
