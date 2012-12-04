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
#include "../../../lattice/access_pattern/StandardPattern.hxx"
#include "../../../lattice/access_pattern/GpuCoulombPatternParity.hxx"
#include "../../../lattice/access_pattern/GpuLandauPatternParity.hxx"
#include "../../../lattice/SiteCoord.hxx"
#include "../../../lattice/SiteIndex.hxx"
#include "../../../lattice/LinkFile.hxx"
#include "../../../util/timer/Chronotimer.h"
#include "../../../lattice/filetypes/FileHeaderOnly.hxx"
#include "../../../lattice/filetypes/FilePlain.hxx"
#include "../../../lattice/filetypes/FileVogt.hxx"
#include "../../../lattice/filetypes/filetype_typedefs.h"
#include "../../GlobalConstants.h"
#include "../program_options/ProgramOptions.hxx"
#include "../program_options/FileIterator.hxx"

// #include "MultiGPU_MPI_LandauGaugeFixingSU3_4D.h"
#include "./MultiGPU_MPI_Communicator.hxx"
#include "./MultiGPU_MPI_AlgorithmOptions.h"


using namespace std;

//TODO where to put these constants?
// const lat_dim_t Ndim = 4;
// const short Nc = 3;
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
// const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

// lattice setup
// const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
// const lat_coord_t sizeTimeslice[Ndim] = {1,Nx,Ny,Nz};



typedef GpuCoulombPatternParity<SiteCoord<Ndim,TIMESLICE_SPLIT>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;

void initNeighbourTable( lat_index_t* nnt )
{
// 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteIndex<4,FULL_SPLIT> s(size);
	SiteIndex<4,FULL_SPLIT> s( HOST_CONSTANTS::SIZE_TIMESLICE);
	s.calculateNeighbourTable( nnt );
}





int main(int argc, char* argv[])
{
	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;
	
	// instantiate object of MPI communicator
	MultiGPU_MPI_Communicator< MultiGPU_MPI_LandauKernelsSU3 > comm(argc,argv);
	
	Chronotimer kernelTimer;
	if( comm.isMaster() ) kernelTimer.reset();
	if( comm.isMaster() ) kernelTimer.start();
	
	// inst. obj. to pass algorithm options to kernel wrappers
	MultiGPU_MPI_AlgorithmOptions algoOptions;
	
	algoOptions.setSaSteps( options.getSaSteps() );
	algoOptions.setSaMin( options.getSaMin() );
	algoOptions.setSaMax( options.getSaMax() );
	algoOptions.setSaMicroupdates( options.getSaMicroupdates() );
	algoOptions.setOrParameter( options.getOrParameter() );
	algoOptions.setSrParameter( options.getSrParameter() );

	Chronotimer allTimer;
	if( comm.isMaster() ) allTimer.reset();

	SiteCoord<4,TIMESLICE_SPLIT> s(HOST_CONSTANTS::SIZE);
	
	// TODO maybe we should choose the filetype at compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfHeaderOnly( options.getReinterpret() );
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfVogt( options.getReinterpret() );
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfPlain( options.getReinterpret() );
	
// 	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfHeaderOnly;
// 	lfHeaderOnly.reinterpret = options.getReinterpret();
// 	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfVogt;
// 	lfVogt.reinterpret = options.getReinterpret();
// 	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfPlain;
// 	lfPlain.reinterpret = options.getReinterpret();
	
	
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

	
// 	TODO cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );
	

	

// 	float totalKernelTime = 0;
// 	long totalStepNumber = 0;
	
	double orTotalKernelTime = 0; // sum up total kernel time for OR
	long orTotalStepnumber = 0;
	double saTotalKernelTime = 0;

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
		//			lfVogt.reinterpret = options.getReinterpret();
					loadOk = lfVogt.load( s, fi.getFilename(), U );
					break;
				case PLAIN:
		//			lfPlain.reinterpret = options.getReinterpret();
					loadOk = lfPlain.load( s, fi.getFilename(), U );
					break;
				case HEADERONLY:
		//			lfHeaderOnly.reinterpret = options.getReinterpret();
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
		//---------------------------------------------------------------------------
		// send (parts of) config. to all processes and transfer to devices
// 		comm.scatterGaugeField( dU, U );
		
		
// 		// calculate and print the gauge quality
// 		if( comm.isMaster() ) printf( "i:\t\tgff:\t\tdA:\n");
// 		comm.generateGaugeQuality( dU, dNnt );
// 		
// 		//TODO divide somewhere else
// 		if( comm.isMaster() ) printf( "-\t\t%1.10f\t\t%e\n", comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
			
			
// 		Chronotimer kernelTimer;
// 		if( comm.isMaster() ) kernelTimer.reset();
// 		if( comm.isMaster() ) kernelTimer.start();
// 		
// 		
// 		// set algorithm = overrelaxation
// 		algoOptions.setAlgorithm( OR );
// 		
// 		// iterate the algorithm
// 		for( int i = 0; i < options.getOrMaxIter(); i++ )
// 		{
// 			for( int evenodd=0; evenodd<2; evenodd++ )
// 			{
// 				comm.apply( dU, dNnt, evenodd, algoOptions );
// 			} 
// 
// 			// calculate and print the gauge quality
// 			if( i % options.getCheckPrecision() == 0 )
// 			{
// 				comm.generateGaugeQuality( dU, dNnt );
// 				
// 				if( comm.isMaster() ) printf( "%d\t\t%1.10f\t\t%e\n", i, comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
// 			}
// 
// 			totalStepNumber++;
// 		}
//  		cudaDeviceSynchronize();
// 		if( comm.isMaster() ) kernelTimer.stop();
// 		if( comm.isMaster() ) cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
// 		totalKernelTime += kernelTimer.getTime();
		
		//---------------------------------------------------------------------------
		
		double bestGff = 0.0;
		for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
		{
			// we copy from host in every gaugecopy step to have a cleaner configuration (concerning numerical errors)
			comm.scatterGaugeField( dU, U );

			// random trafo
			if( !options.isNoRandomTrafo() ) // I'm an optimist! This should be called isRandomTrafo()!
			{
				// set algorithm = random transformation
				algoOptions.setAlgorithm( RT );
				comm.apply( dU, dNnt, 0, algoOptions );
				comm.apply( dU, dNnt, 1, algoOptions );
			}

			// calculate and print the gauge quality
			if( comm.isMaster() ) printf( "i:\t\tgff:\t\tdA:\n");
			comm.generateGaugeQuality( dU, dNnt );
			
			//TODO divide somewhere else
			if( comm.isMaster() ) printf( "-\t\t%1.10f\t\t%e\n", comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
			

			algoOptions.setTemperature( options.getSaMax() );
			algoOptions.setTempStep( (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps() );

			if( comm.isMaster() ) kernelTimer.reset();
			if( comm.isMaster() ) kernelTimer.start();
			
			// simulated annealing
			for( int i = 0; i < options.getSaSteps(); i++ )
			{
				// set algorithm = simulated annealing
				algoOptions.setAlgorithm( SA );
				comm.apply( dU, dNnt, 0, algoOptions );
				comm.apply( dU, dNnt, 1, algoOptions );

				for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
				{
					// set algorithm = micro step
					algoOptions.setAlgorithm( MS );
					comm.apply( dU, dNnt, 0, algoOptions );
					comm.apply( dU, dNnt, 1, algoOptions );
				}


				if( i % options.getCheckPrecision() == 0 )
				{
// 					comm.projectSU3( dU );
					comm.generateGaugeQuality( dU, dNnt );
					if( comm.isMaster() ) printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, algoOptions.getTemperature(), comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );
				}
				algoOptions.decreaseTemperature();
			}
// 			cudaThreadSynchronize();
			if( comm.isMaster() ) kernelTimer.stop();
			if( comm.isMaster() ) saTotalKernelTime += kernelTimer.getTime();

			if( comm.isMaster() ) kernelTimer.reset();
			if( comm.isMaster() ) kernelTimer.start();
			
			// overrelaxation
			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{
				// set algorithm = overrelaxation
				algoOptions.setAlgorithm( OR );
				comm.apply( dU, dNnt, 0, algoOptions );
				comm.apply( dU, dNnt, 1, algoOptions );

				if( i % options.getCheckPrecision() == 0 )
				{
// 					comm.projectSU3( dU );
					comm.generateGaugeQuality( dU, dNnt );
					if( comm.isMaster() ) printf( "%d\t\t%1.10f\t\t%e\n", i, comm.getCurrentGff()/double(s.getLatticeSize())/(double)Ndim/(double)Nc, comm.getCurrentA()/double(s.getLatticeSize())/(double)Nc );;
				}

				if( comm.isMaster() ) orTotalStepnumber++;
			}

// 			cudaThreadSynchronize();
			if( comm.isMaster() ) kernelTimer.stop();
			if( comm.isMaster() ) cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
			if( comm.isMaster() ) orTotalKernelTime += kernelTimer.getTime();




			// reconstruct third line
			comm.projectSU3( dU );

			// check for best copy
			if( comm.getCurrentGff() > bestGff )
			{
				if( comm.isMaster() ) cout << "FOUND BETTER COPY" << endl;
				bestGff = comm.getCurrentGff();

				// send back all timeslices to master
				comm.collectGaugeField( dU, U );
			}
			else
			{
				if( comm.isMaster() ) cout << "NO BETTER COPY" << endl;
			}
			
		} // end for copy
		
		
		
		
		//saving file
// 		if( comm.isMaster() )
// 		{
// 			cout << "saving " << fi.getOutputFilename() << " as " << options.getFType() << endl;
// 		switch( options.getFType() )
// 		{
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

// 	if( comm.isMaster() ) allTimer.stop();
// 	if( comm.isMaster() ) cout << "total time: " << allTimer.getTime() << " s" << endl;
// 	if( comm.isMaster() ) cout << "total kernel time: " << totalKernelTime << " s" << endl;
// 	if( comm.isMaster() ) cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
// 				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

	if( comm.isMaster() )
	{
		long hbFlops = 2176;
		long microFlops = 2118;
		cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
						<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

		long orFlops = 2253;
// 		long orFlops = 2124; // HV 2012-12-03
		cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;
	}
				
	return 0;

}
