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
#include "../lattice/gaugefixing/GaugeFixingStatsV2.hxx"
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

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;



//// boost program options setup
//boost::program_options::variables_map options_vm;
//boost::program_options::options_description options_desc("Allowed options");
//
//// parameters from command line or config file
//int nconf;
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
//int fileNumberformat;
//string configFile;
//bool noRandomTrafo;
//FileType fileType;
//ReinterpretReal reinterpretReal;

// lattice setup
//const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef GpuCoulombPattern<SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> GpuTimeslice;


typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;
typedef Link<GpuTimeslice,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink3;

void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE_TIMESLICE);
	s.calculateNeighbourTable( nnt );
}



__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteCoord<4,FULL_SPLIT> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink3 linkUp( U, s, mu );
		SU3<TLink3> globUp( linkUp );

		globUp.projectSU3();
	}
}

__global__ void projectSU3DP( Real* U )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteCoord<4,FULL_SPLIT> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink3 linkUp( U, s, mu );
		SU3<TLink3> globUp( linkUp );

		Matrix<complex,Nc> locMat;
		SU3<Matrix<complex,Nc> > temp(locMat);

		temp.assignWithoutThirdLine( globUp );

		Matrix<Complex<double>,Nc > doubleMat;

		for( int i = 0; i < 2; i++ )
			for(int j = 0; j < 3; j++ )
			{
				Complex<double> a;
				a.x = (double)temp.get( i, j ).x;
				a.y = (double)temp.get( i, j ).y;
				doubleMat.set( i,j, a);
			}
		double abs_u = 0, abs_v = 0;
		Complex<double> sp(0.,0.);

		// normalize first row
		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			abs_u += doubleMat.get( 0, i ).abs_squared();
		}

		abs_u = sqrt(abs_u);

		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			doubleMat.set( 0, i, doubleMat.get(0,i)/abs_u );
		}

		// orthogonalize second row
		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			sp += doubleMat.get( 1,i ) * doubleMat.get( 0, i ).conj();
		}
		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			doubleMat.set( 1, i, doubleMat.get(1,i) - doubleMat.get( 0, i)*sp );
		}

		// normalize second row
		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			abs_v += doubleMat.get( 1, i ).abs_squared();
		}
		abs_v = sqrt(abs_v);
		for( lat_group_dim_t i = 0; i < 3; i++ )
		{
			doubleMat.set( 1, i, doubleMat.get(1,i)/abs_v );
		}


		for( int i = 0; i < 2; i++ )
			for(int j = 0; j < 3; j++ )
			{
				Complex<Real> a;
				a.x = (Real)doubleMat.get( i, j ).x;
				a.y = (Real)doubleMat.get( i, j ).y;
				temp.set( i,j, a);
			}

		globUp.assignWithoutThirdLine( temp );
	}
}

//__global__ void generateGaugeQuality( Real *U, Real *dGff, Real *dA )
//{
//	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
//	SiteCoord<3,true> s(size);
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	Matrix<Complex<Real>,Nc> locMatSum;
//	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);
//
//	Sum.zero();
//
//	// TODO check if there is a faster way to compute DELTA
//	for( int mu = 1; mu < 4; mu++ )
//	{
//		s.setLatticeIndex( site );
//
//		Matrix<Complex<Real>,Nc> locMat;
//		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);
//
//		TLink3 linkUp( U, s, mu );
//		SU3<TLink3> globUp( linkUp );
//
//		temp.assignWithoutThirdLine( globUp );
//		temp.reconstructThirdLine();
//		Sum += temp;
//
//		s.setNeighbour(mu-1,-1);
//		TLink3 linkDw( U, s, mu );
//		SU3<TLink3> globDw( linkDw );
//		temp.assignWithoutThirdLine( globDw );
//		temp.reconstructThirdLine();
//		Sum -= temp;
//	}
//
//	Sum -= Sum.trace()/Real(3.);
//
//	Matrix<Complex<Real>,Nc> locMatSumHerm;
//	SU3<Matrix<Complex<Real>,Nc> > SumHerm(locMatSumHerm);
//	SumHerm = Sum;
//	SumHerm.hermitian();
//
//	Sum -= SumHerm;
//
//	Real prec = 0;
//	for( int i = 0; i < 3; i++ )
//	{
//		for( int j = 0; j < 3; j++ )
//		{
//			prec += Sum.get(i,j).abs_squared();
//		}
//	}
//
//	dA[site] = prec;
//
//
//	s.setLatticeIndex( site );
//	Real result = 0;
//
//
//	Matrix<Complex<Real>,Nc> locTemp;
//	SU3<Matrix<Complex<Real>,Nc> > temp(locTemp);
//	for( int mu = 1; mu < 4; mu++ )
//	{
//		TLink3 linkUp( U, s, mu );
//		SU3<TLink3> globUp( linkUp );
//		temp.assignWithoutThirdLine( globUp ); // TODO don't load twice
//		temp.reconstructThirdLine();
//		result += temp.trace().x;
//	}
//
//	dGff[site] = result;
//}



Real calculatePolyakovLoopAverage( Real *U )
{
	Matrix<Complex<Real>,3> tempMat;
	SU3<Matrix<Complex<Real>,3> > temp( tempMat );
	Matrix<Complex<Real>,3> temp2Mat;
	SU3<Matrix<Complex<Real>,3> > temp2( temp2Mat );

	SiteCoord<Ndim,FULL_SPLIT> s( HOST_CONSTANTS::SIZE );

	Complex<Real> result(0,0);

	for( s[1] = 0; s[1] < s.size[1]; s[1]++ )
	{
		for( s[2] = 0; s[2] < s.size[2]; s[2]++ )
		{
			for( s[3] = 0; s[3] < s.size[3]; s[3]++ )
			{
				temp.identity();
				temp2.zero();

				for( s[0] = 0; s[0] < s.size[0]; s[0]++ )
				{

					TLink link( U, s, 0 );
					SU3<TLink> globU( link );

					temp2 = temp2 + temp*globU;

					temp = temp2;
					temp2.zero();
				}
				result += temp.trace();
			}
		}
	}

	return sqrt(result.x*result.x+result.y*result.y) / (Real)(s.getLatticeSizeTimeslice()*Nc);
}






int main(int argc, char* argv[])
{
	CoulombKernelsSU3::initCacheConfig();


	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;


	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	printf("\nDevice %d: \"%s\"\n", 0, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);

	Chronotimer allTimer;
	allTimer.reset();

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

//	// device memory for collecting the parts of the gauge fixing functional and divA
//	Real *dGff;
//	cudaMalloc( &dGff, s.getLatticeSizeTimeslice()*sizeof(Real) );
//	Real *dA;
//	cudaMalloc( &dA, s.getLatticeSizeTimeslice()*sizeof(Real) );

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

	allTimer.start();

	GaugeFixingStats<Ndim,Nc,CoulombKernelsSU3,AVERAGE> gaugeStats( dUtUp, HOST_CONSTANTS::SIZE_TIMESLICE );

	double totalKernelTime = 0;
	long totalStepNumber = 0;


	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
	{
//	for( int i = options.getFStartnumber(); i < options.getFStartnumber()+options.getNconf(); i++ )
//	{
//		stringstream filename(stringstream::out);
//		filename << options.getFBasename() << setw( options.getFNumberformat() ) << setfill( '0' ) << i << options.getFEnding();
//		filename << "/home/vogt/configs/STUDIENARBEIT/N32/config_n32t32beta570_sp" << setw( 4 ) << setfill( '0' ) << i << ".vogt";
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

			// copying timeslice t ...
			cudaMemcpy( dUtUp, &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
			// ... and t-1 to device
			cudaMemcpy( dUtDw, &U[tDw*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
			// TODO it is not necessary to copy the (t-1) again for t>0, simply swap pointers on device side...

			// calculate and print the gauge quality
			printf( "i:\t\tgff:\t\tdA:\n");
			gaugeStats.generateGaugeQuality();

			Chronotimer kernelTimer;
			kernelTimer.reset();
			kernelTimer.start();

			float temperature = options.getSaMax();
			float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

			for( int i = 0; i < options.getSaSteps(); i++ )
			{


				CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, temperature, PhiloxWrapper::getNextCounter() );
				CoulombKernelsSU3::saStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, temperature, PhiloxWrapper::getNextCounter() );

				for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
				{
					CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0 );
					CoulombKernelsSU3::microStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1 );
				}


				if( i % options.getOrCheckPrecision() == 0 )
				{
					projectSU3<<<numBlocks*2,32>>>( dUtUp );
					projectSU3<<<numBlocks*2,32>>>( dUtDw );


					gaugeStats.generateGaugeQuality();
					CudaError::getLastError( "generateGaugeQuality error" );
					printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
				}
				temperature -= tempStep;
			}

//			projectSU3DP<<<numBlocks*2,32>>>( dUtUp );
//			projectSU3DP<<<numBlocks*2,32>>>( dUtDw );

			for( int i = 0; i < options.getOrMaxIter(); i++ )
			{

				CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 0, options.getOrParameter() );
				CoulombKernelsSU3::orStep(numBlocks,threadsPerBlock,dUtUp, dUtDw, dNnt, 1, options.getOrParameter() );

//				if( i < 5000 && i % 5 == 0 )
//				{
//					projectSU3<<<numBlocks*2,32>>>( dUtUp );
//					projectSU3<<<numBlocks*2,32>>>( dUtDw );
//				}

				if( i % options.getOrCheckPrecision() == 0 )
				{
					projectSU3<<<numBlocks*2,32>>>( dUtUp );
					projectSU3<<<numBlocks*2,32>>>( dUtDw );
					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < options.getOrPrecision() ) break;
				}

				totalStepNumber++;
			}

			cudaThreadSynchronize();
			kernelTimer.stop();
			cout << "kernel time for timeslice: " << kernelTimer.getTime() << " s"<< endl;
			totalKernelTime += kernelTimer.getTime();
			// copy back TODO: copying back timeslice t is not necessary (only in the end)
			cudaMemcpy( &U[t*timesliceArraySize], dUtUp, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
			cudaMemcpy( &U[tDw*timesliceArraySize], dUtDw, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );

//			exit(-1);
		}


//		filename << ".gf";
//		bool saveOk = lf.save( s, filename.str(), U );
//		if( !saveOk )
//		{
//			cout << "Error while writing." << endl;
//			break;
//		}
//		else
//		{
//			cout << "File written." << endl;
//		}

	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
	cout << "total kernel time: " << totalKernelTime << " s" << endl;
	cout << (double)((long)2205*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9/(double)s.size[0] << " GFlops at "
			<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9/(double)s.size[0] << "GB/s memory throughput." << endl;


}
