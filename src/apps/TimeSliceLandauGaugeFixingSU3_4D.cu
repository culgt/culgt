/*
 * test_gaugefixing.cpp
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt&schroeck
 */

#include <iostream>
#include <math.h>
#include <sstream>
#ifndef OSX
#include "malloc.h"
#endif
#include "../lattice/gaugefixing/GaugeFixingSubgroupStep.hxx"
#include "../lattice/gaugefixing/GaugeFixingStats.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/LinkFile.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/filetype_typedefs.h"
#include "../util/datatype/lattice_typedefs.h"
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>


using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;

#ifdef _X_
const lat_coord_t Nx = _X_;
#else
#error "Define X (the lattice size in x-direction)"
#endif
#ifdef _Y_
const lat_coord_t Ny = _Y_;
#else
const lat_coord_t Ny = _X_;
bool warnY = true; // TODO print the warning
#endif
#ifdef _Z_
const lat_coord_t Nz = _Z_;
#else
const lat_coord_t Nz = _X_;
bool warnZ = true;
#endif
#ifdef _T_
const lat_coord_t Nt = _T_;
#else
#error "Define T (the lattice size in t-direction)"
#endif


// boost program options setup
boost::program_options::variables_map options_vm;
boost::program_options::options_description options_desc("Allowed options");

// parameters from command line or config file
int nconf;
int numbdevices;
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
__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef GpuCoulombPattern<SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,false>,Ndim,Nc> Standard;

typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;

void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.calculateNeighbourTable( nnt );
}


// __device__ inline Real cuFabs( Real a )
// {
// 	return (a>0)?(a):(-a);
// }


__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteCoord<4,true> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink linkUp( U, s, mu );
		SU3<TLink> globUp( linkUp );


		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

		locU.assignWithoutThirdLine(globUp);
		locU.projectSU3withoutThirdRow();


		globUp.assignWithoutThirdLine(locU);

//		globUp.projectSU3withoutThirdRow();
	}
}


__global__ void __launch_bounds__(256,4) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter, int counter=0  )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.nn = nnt;

	const bool updown = threadIdx.x / 128;
	const short mu = (threadIdx.x % 128) / 32;
	const short id = (threadIdx.x % 128) % 32;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );

	Real* U;
	
	//TODO make nicer:
	if( updown==1 )
	{
		if( mu!=0 )
		{
			s.setNeighbour(mu,0);
			U=UtUp;
		}
		else
		{
			U=UtDw;
		}
	}
	else
		U=UtUp;

	Matrix<Complex<Real>,Nc> locMat;
	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
	TLinkIndex link( U, s, mu );
	SU3<TLinkIndex> globU( link );
	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();

	// define the update algorithm
	OrUpdate overrelax( orParameter );
	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, OrUpdate, LANDAU> subgroupStep( &locU, overrelax, id, mu, updown, counter );

	// do the subgroup iteration
	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );

	// project back
// 	globU.projectSU3withoutThirdRow();
	
	// copy link back
	globU.assignWithoutThirdLine(locU);
	//globU=locU; //TODO with or without 3rd line?
}


__global__ void generateGaugeQuality( Real* UtUp, Real* UtDw, lat_index_t* nnt, float *dGff, float *dA )
{
	typedef GpuLandauPattern< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.nn = nnt;
	
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	if( site >= Nx*Ny*Nz ) return; //important in case Nx^3 is not power of 2

	Matrix<Complex<Real>,Nc> locMatSum;
	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);
	Sum.zero();
	float result = 0;


	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > temp(locMat);
		TLinkIndex linkUp( UtUp, s, mu );
		SU3<TLinkIndex> globUp( linkUp );
		
		temp.assignWithoutThirdLine( globUp );
		temp.reconstructThirdLine();
		result += temp.trace().x;
		Sum += temp;
		
		if( mu==0 )
		{
			TLinkIndex linkDw( UtDw, s, mu );
			SU3<TLinkIndex> globDw( linkDw );
			temp.assignWithoutThirdLine( globDw );
			temp.reconstructThirdLine();
		}
		else
		{
			s.setNeighbour(mu,false);
			TLinkIndex linkDw( UtUp, s, mu );
			SU3<TLinkIndex> globDw( linkDw );
			temp.assignWithoutThirdLine( globDw );
			temp.reconstructThirdLine();
		}
		Sum -= temp;
	}
	dGff[site] = result;

	Sum -= Sum.trace()/Real(3.);

	Matrix<Complex<Real>,Nc> locMatSumHerm;
	SU3<Matrix<Complex<Real>,Nc> > SumHerm(locMatSumHerm);
	SumHerm = Sum;
	SumHerm.hermitian();

	Sum -= SumHerm;

	float prec = 0;
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			prec += Sum.get(i,j).abs_squared();
		}
	}

	dA[site] = prec;
}


__global__ void averageGaugeQuality( float* dGff, float* dA )
{
	float gff = 0;
	float A = 0;
	for( int i = 0; i < Nx*Ny*Nz; i++ )
	{
		gff+= dGff[i];
		A += dA[i];
	}
	
	dGff[0] = gff/float(Nx*Ny*Nz)/4./3.;
	dA[0]   = A/float(Nx*Ny*Nz)/3.;
}




// Real calculatePolyakovLoopAverage( Real *U )
// {
// 	Matrix<Complex<Real>,3> tempMat;
// 	SU3<Matrix<Complex<Real>,3> > temp( tempMat );
// 	Matrix<Complex<Real>,3> temp2Mat;
// 	SU3<Matrix<Complex<Real>,3> > temp2( temp2Mat );
// 
// 	SiteCoord<Ndim,true> s( size );
// 
// 	Complex<Real> result(0,0);
// 
// 	for( s[1] = 0; s[1] < s.size[1]; s[1]++ )
// 	{
// 		for( s[2] = 0; s[2] < s.size[2]; s[2]++ )
// 		{
// 			for( s[3] = 0; s[3] < s.size[3]; s[3]++ )
// 			{
// 				temp.identity();
// 				temp2.zero();
// 
// 				for( s[0] = 0; s[0] < s.size[0]; s[0]++ )
// 				{
// 
// 					TLink link( U, s, 0 );
// 					SU3<TLink> globU( link );
// 
// 					temp2 = temp2 + temp*globU;
// 
// 					temp = temp2;
// 					temp2.zero();
// 				}
// 				result += temp.trace();
// 			}
// 		}
// 	}
// 
// 	return sqrt(result.x*result.x+result.y*result.y) / (Real)(s.getLatticeSizeTimeslice()*Nc);
// }


int main(int argc, char* argv[])
{
	// read parameters (command line or given config file)
	options_desc.add_options()
		("help", "produce help message")
		("numbdevices", boost::program_options::value<int>(&numbdevices)->default_value(1), "how many CUDA devices to use")
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


	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	deviceCount = ( deviceCount<numbdevices ? deviceCount : numbdevices );
// 	int device;
// 	for (device = 0; device < deviceCount; ++device) 
// 	{
// 		cudaDeviceProp deviceProp;
    
// 		cudaGetDeviceProperties(&deviceProp, device);
    
// 		printf("Device %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
// 	}
	
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	printf("\nDevice %d: \"%s\"\n", 0, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);


	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,true> s(size);

	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,true> > lfHeaderOnly;
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,true> > lfVogt;
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,true> > lfPlain;


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for all timeslices
	Real* dU[Nt];
	for( int t=0; t<Nt; t++ )
		cudaMalloc( &dU[t], timesliceArraySize*sizeof(Real) );

	// device memory for collecting the parts of the gauge fixing functional and divA
	float *dGff;
	cudaMalloc( &dGff, s.getLatticeSizeTimeslice()*sizeof(float) );
	float *dA;
	cudaMalloc( &dA, s.getLatticeSizeTimeslice()*sizeof(float) );
	float gff[2], A[2];

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt;
	cudaMalloc( &dNnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	// copy neighbour table to device
	cudaMemcpy( dNnt, nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t), cudaMemcpyHostToDevice );

	int threadsPerBlock = 32*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSizeTimeslice()/2/32; // // half of the lattice sites (a parity) are updated in a kernel call

	allTimer.start();

	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );

// 	lat_coord_t *pointerToSize;
// 	cudaGetSymbolAddress( (void**)&pointerToSize, "dSize" );

	float totalKernelTime = 0;
	long totalStepNumber = 0;


	for( int i = fileStartnumber; i < fileStartnumber+nconf; i++ )
	{
		stringstream filename(stringstream::out);
		filename << fileBasename << setw( fileNumberformat ) << setfill( '0' ) << i << fileEnding;
//		filename << "/home/vogt/configs/STUDIENARBEIT/N32/config_n32t32beta570_sp" << setw( 4 ) << setfill( '0' ) << i << ".vogt";
		cout << "loading " << filename.str() << " as " << fileType << endl;
		cout << filename.str() << endl;
		bool loadOk;

		switch( fileType )
		{
		case VOGT:
			loadOk = lfVogt.load( s, filename.str(), U );
			break;
		case PLAIN:
			loadOk = lfPlain.load( s, filename.str(), U );
			break;
		case HEADERONLY:
			loadOk = lfHeaderOnly.load( s, filename.str(), U );
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
		
// 		Real polBefore = calculatePolyakovLoopAverage( U );

		// copying all timeslices to devices
		for( int t=0; t<Nt; t++ )
			cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );

		// calculate and print the gauge quality
		printf( "i:\t\tgff:\t\tdA:\n");
		gff[1]=0.0; 
		A[1]  =0.0;
		for( int t=0; t<Nt; t++ )
		{
// 			cout << t << ' ' << flush;
			int tDw = (t > 0)?(t-1):(s.size[0]-1);
			generateGaugeQuality<<<Nx*Nx*Nx/32,32>>>( dU[t], dU[tDw], dNnt, dGff, dA );
			averageGaugeQuality<<<1,1>>>( dGff, dA );
			cudaMemcpy( gff, dGff, sizeof(float), cudaMemcpyDeviceToHost );
			cudaMemcpy( A, dA,     sizeof(float), cudaMemcpyDeviceToHost );
			gff[1]+=gff[0];
			A[1]  +=A[0];
		}
		printf( "\n-\t\t%1.10f\t\t%e\n", gff[1]/(float)Nt, A[1]/(float)Nt );
		
		Chronotimer kernelTimer;
		kernelTimer.reset();
		kernelTimer.start();
		
		for( int j = 0; j < orMaxIter; j++ )
		{
			for( int t=0; t<Nt; t++ )
			{
// 				cout << t << endl;
				int tDw = (t > 0)?(t-1):(s.size[0]-1);
				orStep<<<numBlocks,threadsPerBlock>>>(dU[t], dU[tDw], dNnt, 0, orParameter );
				orStep<<<numBlocks,threadsPerBlock>>>(dU[t], dU[tDw], dNnt, 1, orParameter );
			}

			// calculate and print the gauge quality
			if( j % orCheckPrec == 0 )
			{
				gff[1]=0.0; 
				A[1]  =0.0;
				for( int t=0; t<Nt; t++ )
				{
// 					cout << t << ' ' << flush;
					int tDw = (t > 0)?(t-1):(s.size[0]-1);
					generateGaugeQuality<<<Nx*Nx*Nx/32,32>>>( dU[t], dU[tDw], dNnt, dGff, dA );
					averageGaugeQuality<<<1,1>>>( dGff, dA );
					cudaMemcpy( gff, dGff, sizeof(float), cudaMemcpyDeviceToHost );
					cudaMemcpy( A, dA,     sizeof(float), cudaMemcpyDeviceToHost );
					gff[1]+=gff[0];
					A[1]  +=A[0];
				}
				printf( "%d\t\t%1.10f\t\t%e\n", j, gff[1]/(float)Nt, A[1]/(float)Nt );
				
				if( A[1] < orPrecision ) break;
			}

			totalStepNumber++;
		}
		cudaThreadSynchronize();
		kernelTimer.stop();
		cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		
// 		cout << "Polyakov loop: " << polBefore << " - " << calculatePolyakovLoopAverage( U ) << endl;

		// copy back all timeslices
		for( int t=0; t<Nt; t++ )
			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );

	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
	cout << "total kernel time: " << totalKernelTime << " s" << endl;
	cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;



}
