/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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
#include "../../../lattice/access_pattern/GpuPatternTimesliceParityPriority.hxx"
#include "../../../lattice/access_pattern/GpuPatternParityPriority.hxx"
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



typedef GpuPatternTimesliceParityPriority<SiteCoord<Ndim,TIMESLICE_SPLIT>,Ndim,Nc> Gpu;
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
	algoOptions.setSeed( options.getSeed() + comm.getRank() );

	Chronotimer allTimer;
	if( comm.isMaster() ) allTimer.reset();

	SiteCoord<4,TIMESLICE_SPLIT> s(HOST_CONSTANTS::SIZE);
	
	// TODO maybe we should choose the filetype at compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfHeaderOnly( options.getReinterpret() );
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfVogt( options.getReinterpret() );
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,TIMESLICE_SPLIT> > lfPlain( options.getReinterpret() );
	
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
		if( comm.isMaster() && !options.isSetHot() )
		{
			cout << "loading " << fi.getFilename() << " as " << options.getFType() << endl;
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

		// don't read gauge field, set hot:
		if( options.isSetHot() ) comm.setHot( dU, algoOptions );


		
		double bestGff = 0.0;
		for( int copy = 0; copy < options.getGaugeCopies(); copy++ )
		{
			// we copy from host in every gaugecopy step to have a cleaner configuration (concerning numerical errors)
			if( !options.isSetHot() ) comm.scatterGaugeField( dU, U );

			// random trafo
			if( options.isRandomTrafo() )
			{
				// set algorithm = random transformation
				algoOptions.setAlgorithm( RT );
				comm.apply( dU, dNnt, 0, algoOptions );
				comm.apply( dU, dNnt, 1, algoOptions );
			}

			// calculate and print the gauge quality
			if( comm.isMaster() ) printf( "i:\t\tgff:\t\tdA:\n");
			comm.generateGaugeQuality( dU, dNnt );
			
			//print the gauge quality
			if( comm.isMaster() ) printf( "-\t\t%1.10f\t\t%e\n", comm.getCurrentGff(), comm.getCurrentA() );
			

			algoOptions.setTemperature( options.getSaMax() );
			algoOptions.setTempStep( (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps() );

			if( comm.isMaster() ) kernelTimer.reset();
			if( comm.isMaster() ) kernelTimer.start();
			
			
			// SIMULATED ANNEALING
			if( options.getSaSteps()>0  && comm.isMaster() ) 
				printf( "SIMULATED ANNEALING\n" );
			
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

				if( i % options.getReproject() == 0 )
				{
					comm.projectSU3( dU );
				}

				if( i % options.getCheckPrecision() == 0 )
				{
					comm.generateGaugeQuality( dU, dNnt );
					if( comm.isMaster() ) printf( "%d\t%f\t\t%1.10f\t\t%e\n", i, algoOptions.getTemperature(), comm.getCurrentGff(), comm.getCurrentA() );
					if( comm.getCurrentA() < options.getPrecision() ) break;
				}
				algoOptions.decreaseTemperature();
			}
// 			cudaThreadSynchronize();
			if( comm.isMaster() ) 
			{
				kernelTimer.stop();
				cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
				saTotalKernelTime += kernelTimer.getTime();

				kernelTimer.reset();
				kernelTimer.start();
			}
			
			
			// OVERRELAXATION
			if( options.getOrMaxIter()>0  && comm.isMaster() ) 
				printf( "OVERRELAXATION\n" );
			
			// set algorithm = overrelaxation
			algoOptions.setAlgorithm( OR );
			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{
				comm.apply( dU, dNnt, 0, algoOptions );
				comm.apply( dU, dNnt, 1, algoOptions );

				if( i % options.getReproject() == 0 )
				{
					comm.projectSU3( dU );
				}

				if( i % options.getCheckPrecision() == 0 )
				{
					comm.generateGaugeQuality( dU, dNnt );
					if( comm.isMaster() ) printf( "%d\t\t%1.10f\t\t%e\n", i, comm.getCurrentGff(), comm.getCurrentA() );
					if( comm.getCurrentA() < options.getPrecision() ) break;
				}

				if( comm.isMaster() ) orTotalStepnumber++;
			}

// 			cudaThreadSynchronize();
			if( comm.isMaster() ) 
			{
				kernelTimer.stop();
				cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
				orTotalKernelTime += kernelTimer.getTime();
			}




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
		if( comm.isMaster() && !options.isSetHot() )
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
	} // end fileIterator


	if( comm.isMaster() )
	{
		long hbFlops = 2252+86;
		long microFlops = 2252+14;
		cout << "Simulated Annealing (HB+Micro): " << (double)((long)(hbFlops+microFlops*options.getSaMicroupdates())*(long)s.getLatticeSize()*(long)options.getSaSteps()*(long)options.getGaugeCopies())/saTotalKernelTime/1.0e9 << " GFlops at "
						<< (double)((long)192*(long)s.getLatticeSize()*options.getSaSteps()*(options.getSaMicroupdates()+1)*(long)sizeof(Real))/saTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;


		long orFlops = 2252+22;
	//	long orFlops = 2124; // HV 2012-12-03
		cout << "Overrelaxation: " << (double)((long)orFlops*(long)s.getLatticeSize()*(long)orTotalStepnumber)/orTotalKernelTime/1.0e9 << " GFlops at "
					<< (double)((long)192*(long)s.getLatticeSize()*(long)(orTotalStepnumber)*(long)sizeof(Real))/orTotalKernelTime/1.0e9 << "GB/s memory throughput." << endl;
	}
				
	return 0;

}
