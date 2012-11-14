/*
 * test_gaugefixing.cpp
 *
 *  Created on: Nov. 12, 2012
 *      Author: vogt&schroeck
 */

#include <iostream>
#include <math.h>
#include <sstream>
#ifndef OSX
#include "malloc.h"
#endif
// #include "../../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
// #include "../../lattice/gaugefixing/GaugeFixingStats.hxx"
// #include "../../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../../lattice/access_pattern/StandardPattern.hxx"
#include "../../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../../lattice/access_pattern/GpuLandauPatternParity.hxx"
#include "../../lattice/SiteCoord.hxx"
#include "../../lattice/SiteIndex.hxx"
#include "../../lattice/Link.hxx"
#include "../../lattice/SU3.hxx"
#include "../../lattice/Matrix.hxx"
// #include "../../lattice/LinkFile.hxx"
// #include "../../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../../util/timer/Chronotimer.h"
#include "../../lattice/filetypes/FileHeaderOnly.hxx"
#include "../../lattice/filetypes/FilePlain.hxx"
#include "../../lattice/filetypes/FileVogt.hxx"
#include "../../lattice/filetypes/filetype_typedefs.h"
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include "../../lattice/gaugefixing/GlobalConstants.hxx"
// #include "../../util/rng/PhiloxWrapper.hxx"

#include <cuda_runtime.h>
#include "MultiGPU_MPI_LandauGaugeFixingSU3_4D.h"

// MPI include
#include <mpi.h>

// User include
// #include "simpleMPI.h"

// Error handling macros
#define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
        cerr << "MPI error calling \""#call"\"\n"; \
        my_abort(-1); }

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;


// boost program options setup
boost::program_options::variables_map options_vm;
boost::program_options::options_description options_desc("Allowed options");

// parameters from command line or config file
int nconf;
int numbDevices;
long seed; // TODO check datatype
int orMaxIter;
int orCheckPrec;
float orParameter;
float orPrecision;
int saSteps;
float saMin;
float saMax;
int gaugeCopies;
string fileEnding;
string postFixLabel;
string fileBasename;
int fileStartnumber;
int fileStepsize;
int fileNumberformat;
string configFile;
bool noRandomTrafo;
FileType fileType;
ReinterpretReal reinterpretReal;


// lattice setup
const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
const lat_coord_t sizeTimeslice[Ndim] = {1,Nx,Ny,Nz};
#ifdef CUDA
__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
#endif
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef GpuCoulombPattern<SiteCoord<Ndim,TIMESLICE_SPLIT>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
// typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;

void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,TIMESLICE_SPLIT> s(size);
	s.calculateNeighbourTable( nnt );
}



int main(int argc, char* argv[])
{
	// initialize MPI communication
	int nprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  printf("Process %d on %s out of %d alive.\n", rank, processor_name, nprocs);

	int lRank = ( rank - 1 + nprocs ) % nprocs;
	int rRank = ( rank + 1 ) % nprocs;
	bool isMaster = ( rank == 0 ? true : false );
	
	
	
	
	// copy test
// 	float source, sink;
// 	
// 	if( rank==0 ) source=1.2345;
// 	
// // 	double *dS;
// 	
// // 	cudaMalloc( (void**) &dS, sizeof(double) );
// // 	cudaMemcpy( dS, &source, sizeof(double), cudaMemcpyHostToDevice );
// // 	cudaMemcpy( &sink, dS, sizeof(double), cudaMemcpyDeviceToHost );
// 	
// 	MPI_Request request1, request2;
// 	MPI_Status  status;
// 	
// 	MPI_Irecv( &sink,   1, MPI_FLOAT, lRank, 123, MPI_COMM_WORLD, &request2);
// 	MPI_Isend( &source, 1, MPI_FLOAT, rRank, 123, MPI_COMM_WORLD, &request1);
// 	
// 	MPI_Wait( &request1, &status );
// 	MPI_Wait( &request2, &status );
// 	
// 	
// 	cout << sink << endl;
// 	
	
	
	
	
	

//   MPI_Finalize();
// 	return 0;
	
	// read parameters (command line or given config file)
	options_desc.add_options()
		("help", "produce help message")
		("numbdevices", boost::program_options::value<int>(&numbDevices)->default_value(1), "how many CUDA devices to use")
		("nconf,m", boost::program_options::value<int>(&nconf)->default_value(1), "how many files to gaugefix")
		("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(1000), "Max. number of OR iterations")
		("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")
		("sasteps", boost::program_options::value<int>(&saSteps)->default_value(1000), "number of SA steps")
		("samin", boost::program_options::value<float>(&saMin)->default_value(.01), "min. SA temperature")
		("samax", boost::program_options::value<float>(&saMax)->default_value(.4), "max. SA temperature")
		("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "OR parameter")
		("orprecision", boost::program_options::value<float>(&orPrecision)->default_value(1E-7), "OR precision (dmuAmu)")
		("orcheckprecision", boost::program_options::value<int>(&orCheckPrec)->default_value(100), "how often to check the gauge precision")
		("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "Number of gauge copies")
		("ending", boost::program_options::value<string>(&fileEnding)->default_value(".vogt"), "file ending to append to basename")
		("postfixlabel", boost::program_options::value<string>(&postFixLabel)->default_value("_Landau"), "label to append to basename after fixing the gauge and before storing it")
		("basename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
		("startnumber", boost::program_options::value<int>(&fileStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
		("stepsize", boost::program_options::value<int>(&fileStepsize)->default_value(1), "file numbering startnumber, startnumber+stepsize,...")
		("numberformat", boost::program_options::value<int>(&fileNumberformat)->default_value(1), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
		("filetype", boost::program_options::value<FileType>(&fileType), "type of configuration (PLAIN, HEADERONLY, VOGT)")
		("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")
		("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpretReal)->default_value(STANDARD), "reinterpret Real datatype (STANDARD = do nothing, FLOAT = convert input as float and cast to Real, DOUBLE = ...)")
		("norandomtrafo", boost::program_options::value<bool>(&noRandomTrafo)->default_value(false), "no random gauge trafo" )
		;

	boost::program_options::positional_options_description options_p;
	options_p.add("config-file", -1);

	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
			options(options_desc).positional(options_p).run(), options_vm);
	boost::program_options::notify(options_vm);

	ifstream cfg( configFile.c_str() );
	boost::program_options::store(boost::program_options::parse_config_file( cfg, options_desc), options_vm);
	boost::program_options::notify(options_vm);

	if (options_vm.count("help")) {
		cout << "Usage: " << argv[0] << " [options] [config-file]" << endl;
		cout << options_desc << "\n";
		return 1;
	}


	//TODO get rid of:
	int deviceCount;
	deviceCount = ( deviceCount<numbDevices ? deviceCount : numbDevices );
	int device;
	int theDevice[Nt];
	int minT[deviceCount];
	int maxT[deviceCount];
	int midT[deviceCount];

	initDevice( rank );
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	// distribute the timeslices over the threads: get tmin and tmax
	int tmin = rank*Nt/nprocs;
	int tmax = (rank+1)*Nt/nprocs;
	
// 	printf("Process %d: tmin = %d, tmax = %d\n", rank, tmin, tmax);
// 	MPI_Finalize();
// 	return 0;
	
			
	// how many timeslices does each device have?
	int numbTimeSlices=0;
	for( int t=tmin; t<tmax; t++ )
		numbTimeSlices++;
	
	printf("Process %d: numbTimeSlices = %d\n", rank, numbTimeSlices);
	MPI_Barrier(MPI_COMM_WORLD);
    

	
// 	// cudaStreams
// 	cudaStream_t streamStd[deviceCount];
// 	cudaStream_t streamCpy_1[deviceCount];
// 	cudaStream_t streamCpy_2[deviceCount];
// 	for( device=0; device<deviceCount; device++ )
// 	{
// 		cudaSetDevice(device);
// 		cudaStreamCreate(&streamStd[device]);
// 		cudaStreamCreate(&streamCpy_1[device]);
// 		cudaStreamCreate(&streamCpy_2[device]);
// 	}
// 	
// 	// cudaEvents
// 	cudaEvent_t eventCpyD2H[deviceCount];
// 	cudaEvent_t eventCpyH2D[deviceCount];
// 	cudaEvent_t hitFirstTimeslice[deviceCount];
// 	for( device=0; device<deviceCount; device++ )
// 	{
// 		cudaSetDevice(device);
// 		cudaEventCreate(&eventCpyD2H[device]);
// 		cudaEventCreate(&eventCpyH2D[device]);
// 		cudaEventCreate(&hitFirstTimeslice[device]);
// 	}

	Chronotimer allTimer;
	if( isMaster ) allTimer.reset();

	SiteCoord<4,TIMESLICE_SPLIT> s(size);
	
	// TODO maybe we should choose the filetype on compile time
// 	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,true> > lfHeaderOnly;
// 	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,true> > lfVogt;
// 	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,true> > lfPlain;


	// allocate Memory
	// host memory for configuration
// 	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for all timeslices
	Real* dU[Nt];
	for( int t=tmin; t<tmax; t++ )
	{
		cudaMalloc( (void**) &dU[t], timesliceArraySize*sizeof(Real) );
	}
	
	// halo exchange size (only mu=0 and row 1 and 2 of SU(3))
	size_t haloSize = timesliceArraySize*sizeof(Real);
		
	// page-locked host memory for halo timeslices (two per thread)
	Real* haloOut[nprocs];
	Real* haloIn[nprocs];
// 	haloOut[rank] = (Real*)malloc( haloSize );
// 	haloIn[rank] = (Real*)malloc( haloSize );
 	cudaHostAlloc( (void**) &haloIn[rank],  haloSize, 0 );
	cudaHostAlloc( (void**) &haloOut[rank], haloSize, 0 );

	
	// device memory for halo timeslice (one per device)
	Real* dHalo[nprocs];
	cudaMalloc( (void**) &dHalo[rank], haloSize );
	
	
	// device memory for collecting the parts of the gauge fixing functional and divA
	double *dGff[nprocs];
	double *dA[nprocs];
	cudaMalloc( (void**) &dGff[rank], s.getLatticeSizeTimeslice()*sizeof(double) );
	cudaMalloc( (void**) &dA[rank], s.getLatticeSizeTimeslice()*sizeof(double) );


	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt[nprocs];
	cudaMalloc( (void**) &dNnt[rank], s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );


	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	
	// copy neighbour table to device
	cudaMemcpy( dNnt[rank], nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t), cudaMemcpyHostToDevice );
	

	if( isMaster ) allTimer.start();

// 	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );

// 	lat_coord_t *pointerToSize;
// 	cudaGetSymbolAddress( (void**)&pointerToSize, "dSize" );

	float totalKernelTime = 0;
	long totalStepNumber = 0;


	for( int i = fileStartnumber; i < fileStartnumber+nconf; i++ )
	{
// 		stringstream filename(stringstream::out);
// 		filename << fileBasename << setw( fileNumberformat ) << setfill( '0' ) << i << fileEnding;
// //		filename << "/home/vogt/configs/STUDIENARBEIT/N32/config_n32t32beta570_sp" << setw( 4 ) << setfill( '0' ) << i << ".vogt";
// 		cout << "loading " << filename.str() << " as " << fileType << endl;
// 		cout << filename.str() << endl;
// 		bool loadOk;
// 
// 		switch( fileType )
// 		{
// 		case VOGT:
// 			loadOk = lfVogt.load( s, filename.str(), U );
// 			break;
// 		case PLAIN:
// 			loadOk = lfPlain.load( s, filename.str(), U );
// 			break;
// 		case HEADERONLY:
// 			loadOk = lfHeaderOnly.load( s, filename.str(), U );
// 			break;
// 		default:
// 			cout << "Filetype not set to a known value. Exiting";
// 			exit(1);
// 		}
// 
// 		if( !loadOk )
// 		{
// 			cout << "Error while loading. Trying next file." << endl;
// 			break;
// 		}
// 		else
// 		{
// 			cout << "File loaded." << endl;
// 		}
// 		// copying all timeslices to devices
// 		for( int t=0; t<Nt; t++ )
// 		{
// 			cudaSetDevice(rank);
// 			cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
// 		}		
		
		
		// don't load file, fill with random numbers - then projectSU3
		if( isMaster ) cout << "\nWe fill dU with random SU(3) matrices.\n" << endl;
		int unsigned counter = 0;//time(NULL);
		
		for( int t=tmin; t<tmax; t++ )
		{
 			_set_hot( dU[t], counter ); 
			counter++;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		// calculate and print the gauge quality
		if( isMaster ) printf( "i:\t\tgff:\t\tdA:\n");
		// host memory (two fields to collect the results)
		double gff[2], A[2];
		double gffRes, ARes;
		gff[1]=0.0; 
		A[1]  =0.0;
	
		//exchange halos
		MPI_Request request1, request2;
		MPI_Status  status;
		
		cudaMemcpyAsync( haloOut[rank], dU[tmax-1], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Irecv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2);		
		MPI_Isend( haloOut[rank], timesliceArraySize, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1);
		MPI_Wait( &request1, &status );
		MPI_Wait( &request2, &status );
		
		cudaMemcpyAsync( dHalo[rank], haloIn[rank], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
		
		// call kernel wrapper for all t
		for( int t=tmin; t<tmax; t++ )
		{
			if( t == tmin )
				_generateGaugeQualityPerSite( dU[tmin], dHalo[rank], dNnt[rank], dGff[rank], dA[rank] );
			else
				_generateGaugeQualityPerSite( dU[t], dU[t-1], dNnt[rank], dGff[rank], dA[rank] );
			
			_averageGaugeQuality( dGff[rank], dA[rank] );
			cudaMemcpyAsync( &gff[0], dGff[rank], sizeof(double), cudaMemcpyDeviceToHost );
			cudaMemcpyAsync( &A[0], dA[rank],     sizeof(double), cudaMemcpyDeviceToHost );
			gff[1]+=gff[0];
			A[1]  +=A[0];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Reduce( &gff[1], &gffRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( &A[1],   &ARes,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
// 		gffRes = gff[1];
// 		ARes = A[1];
		
		gffRes /= (double)Nt;
		ARes   /= (double)Nt;
		
		if( isMaster ) printf( "-\t\t%1.10f\t\t%e\n", gffRes, ARes );
		
// 		MPI_Finalize();
// 		return 0;
			
		Chronotimer kernelTimer;
		if( isMaster ) kernelTimer.reset();
		if( isMaster ) kernelTimer.start();
		
		for( int j = 0; j < orMaxIter; j++ )
		{
			for( int parity=0; parity<2; parity++ )
			{
				int p_offset = parity?timesliceArraySize/2:0;
					
				// ------------------ ASYNCHRONOUS HALO EXCHANGE --------------------- //
				// copy the halo: device -> host
				cudaMemcpyAsync( haloOut[rank]+p_offset, dU[tmax-1]+p_offset, haloSize/12, cudaMemcpyDeviceToHost );
				MPI_Barrier(MPI_COMM_WORLD);
				// MPI send
				MPI_Irecv( haloIn[rank]+p_offset,  timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2);		
				MPI_Isend( haloOut[rank]+p_offset, timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1);
				MPI_Wait( &request1, &status );
				MPI_Wait( &request2, &status );
				// copy the halo: host -> device
				cudaMemcpyAsync( dHalo[rank]+p_offset, haloIn[rank]+p_offset, haloSize/12, cudaMemcpyHostToDevice );
						
				
				// ------------- HIT TIMESLICES tmin+1,...,(tmax-tmin)/2-1 ----------- //
				for( int t=tmin+1; t<(tmax-tmin)/2; t++ )
				{
// 							cudaEventSynchronize( hitFirstTimeslice[rank] );
					_orStep( dU[t], dU[t-1], dNnt[rank], parity, orParameter );
				}
				
				// -------------------- HIT TIMESLICES tmin ------------------------- //
				// now call kernel wrapper with dU[t-1] replaced by dHalo[rank]
				MPI_Barrier(MPI_COMM_WORLD);
				_orStep( dU[tmin], dHalo[rank], dNnt[rank], parity, orParameter );
				MPI_Barrier(MPI_COMM_WORLD);
				
				// ------------------ ASYNCHRONOUS HALO EXCHANGE --------------------- //
				// copy halo back: device -> host
				cudaMemcpyAsync( haloOut[rank]+p_offset, dHalo[rank]+p_offset, haloSize/12, cudaMemcpyDeviceToHost );
// 			 cudaEventRecord( eventCpyD2H[rank], streamCpy_2[rank] );
				MPI_Barrier(MPI_COMM_WORLD);
				// MPI send
				MPI_Irecv( haloIn[rank]+p_offset,  timesliceArraySize/12, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request2);		
				MPI_Isend( haloOut[rank]+p_offset, timesliceArraySize/12, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request1);
				MPI_Wait( &request1, &status );
				MPI_Wait( &request2, &status );
				// copy halo back: host -> device
// 			 cudaEventSynchronize( eventCpyD2H[rank] );
				cudaMemcpyAsync( dU[tmax-1]+p_offset, haloIn[rank]+p_offset, haloSize/12, cudaMemcpyHostToDevice );
					
				
				// ------------- HIT TIMESLICES (tmax-tmin)/2,...,tmax-1 ----------- //
				for( int t=(tmax-tmin)/2; t<tmax; t++ )
				{
// 				 cudaEventSynchronize( hitFirstTimeslice[rank] );
					_orStep( dU[t], dU[t-1], dNnt[rank], parity, orParameter );
				}
				MPI_Barrier(MPI_COMM_WORLD);
// 				cudaDeviceSynchronize();
			}

			// calculate and print the gauge quality
			if( j % orCheckPrec == 0 )
			{
				gff[1]=0.0; 
				A[1]  =0.0;
			
				cudaMemcpy( haloOut[rank], dU[tmax-1], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
				MPI_Barrier(MPI_COMM_WORLD);
				
				MPI_Irecv( haloIn[rank],  timesliceArraySize, MPI_FLOAT, lRank, 0, MPI_COMM_WORLD, &request2);		
				MPI_Isend( haloOut[rank], timesliceArraySize, MPI_FLOAT, rRank, 0, MPI_COMM_WORLD, &request1);
				MPI_Wait( &request1, &status );
				MPI_Wait( &request2, &status );
				
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
				MPI_Barrier(MPI_COMM_WORLD);
				
				MPI_Reduce( &gff[1], &gffRes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce( &A[1],   &ARes,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		// 		gffRes = gff[1];
		// 		ARes = A[1];
				
				gffRes /= (double)Nt;
				ARes   /= (double)Nt;
				
				if( isMaster ) printf( "%d\t\t%1.10f\t\t%e\n", j, gffRes, ARes );
			}

			totalStepNumber++;
		}
// 		cudaDeviceSynchronize();
		if( isMaster ) kernelTimer.stop();
		if( isMaster ) cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		
// 		cout << "Polyakov loop: " << polBefore << " - " << calculatePolyakovLoopAverage( U ) << endl;

		// copy back all timeslices
// 		for( int t=0; t<Nt; t++ )
// 			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );

	}

	if( isMaster ) allTimer.stop();
	if( isMaster ) cout << "total time: " << allTimer.getTime() << " s" << endl;
	if( isMaster ) cout << "total kernel time: " << totalKernelTime << " s" << endl;
	if( isMaster ) cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


	MPI_Finalize();
	return 0;

}
