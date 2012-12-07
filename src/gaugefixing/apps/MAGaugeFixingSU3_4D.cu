/**
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
 */

#include <iostream>
#include <math.h>
#include <sstream>
#ifndef OSX
#include "malloc.h"
#endif
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

	cudaSetDevice(options.getDeviceNumber());
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, options.getDeviceNumber() );

	printf("\nDevice %d: \"%s\"\n", options.getDeviceNumber(), deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);


	// SiteCoord is faster than SiteIndex when loading files
	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


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
	SiteIndex<4,FULL_SPLIT> sTemp(HOST_CONSTANTS::SIZE);
	sTemp.calculateNeighbourTable( nn );
	// copy neighbour table to device
	cudaMemcpy( dNn, nn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );



	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly( options.getReinterpret() );
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt( options.getReinterpret() );
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain( options.getReinterpret() );


	int threadsPerBlock = NSB*8; // NSB sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSize()/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

	
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
			CommonKernelsSU3::setHot( s.getLatticeSize()/32,32, dU, HOST_CONSTANTS::getPtrToDeviceSize(), options.getSeed(), PhiloxWrapper::getNextCounter() );
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


			// SIMULATED ANNEALING
			if( options.getSaSteps() > 0 ) printf( "SIMULATED ANNEALING\n" );
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

				if( i % options.getReproject() == 0 )
				{
					CommonKernelsSU3::projectSU3( s.getLatticeSize()/32,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
				}

				if( i % options.getCheckPrecision() == 0 )
				{
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
				}
				temperature -= tempStep;
			}
			cudaDeviceSynchronize();
			kernelTimer.stop();
			saTotalKernelTime += kernelTimer.getTime();


			// STOCHASTIC RELAXATION
			if( options.getOrMaxIter() > 0 ) printf( "STOCHASTIC RELAXATION\n" );
			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getSrMaxIter(); i++ )
			{

				MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getSrParameter(), options.getSeed(), PhiloxWrapper::getNextCounter() );
				MAGKernelsSU3::srStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getSrParameter(), options.getSeed(), PhiloxWrapper::getNextCounter() );

				if( i % options.getReproject() == 0 )
				{
					CommonKernelsSU3::projectSU3( s.getLatticeSize()/32,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
				}

				if( i % options.getCheckPrecision() == 0 )
				{
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < options.getPrecision() ) break;
				}
				srTotalStepnumber++;
			}
			cudaDeviceSynchronize();
			kernelTimer.stop();
			srTotalKernelTime += kernelTimer.getTime();


			// OVERRELAXATION
			if( options.getOrMaxIter() > 0 ) printf( "OVERRELAXATION\n" );
			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{

				MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getOrParameter() );
				MAGKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getOrParameter() );

				if( i % options.getReproject() == 0 )
				{
					CommonKernelsSU3::projectSU3( s.getLatticeSize()/32,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );
				}

				if( i % options.getCheckPrecision() == 0 )
				{
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
			CommonKernelsSU3::projectSU3( s.getLatticeSize()/32,32, dU, HOST_CONSTANTS::getPtrToDeviceSize() );

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
