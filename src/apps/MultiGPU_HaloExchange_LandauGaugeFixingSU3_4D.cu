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
#include "../lattice/access_pattern/GpuLandauPatternParity.hxx"
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
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include "../lattice/gaugefixing/GlobalConstants.hxx"
#include "../util/rng/PhiloxWrapper.hxx"


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
__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

// typedef GpuCoulombPattern<SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,false>,Ndim,Nc> Standard;
// typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;

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


// __global__ void projectSU3( Real* U )
// {
// 	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
// 	SiteCoord<4,true> s(size);
// 	int site = blockIdx.x * blockDim.x + threadIdx.x;
// 
// 	s.setLatticeIndex( site );
// 
// 	for( int mu = 0; mu < 4; mu++ )
// 	{
// 		TLink linkUp( U, s, mu );
// 		SU3<TLink> globUp( linkUp );
// 
// 
// 		Matrix<Complex<Real>,Nc> locMat;
// 		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
// 		
// 		locU.assignWithoutThirdLine(globUp);
// 		locU.projectSU3withoutThirdRow();
// 		globUp.assignWithoutThirdLine(locU);
// 	}
// }


__global__ void __launch_bounds__(256,4) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter, int counter=0  )
{
	typedef GpuLandauPatternParity< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.nn = nnt;

	const bool updown = threadIdx.x / NSB4;
	const short mu = (threadIdx.x % NSB4) / NSB;
	const short id = (threadIdx.x % NSB4) % NSB;

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


__global__ void generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, double *dGff, double *dA )
{
	typedef GpuLandauPatternParity< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	s.nn = nnt;
	
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	if( site >= Nx*Ny*Nz ) return; //important in case Nx^3 is not power of 2

	Matrix<Complex<Real>,Nc> locMatSum;
	SU3<Matrix<Complex<Real>,Nc> > Sum(locMatSum);
	Sum.zero();
	double result = 0;


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

	double prec = 0;
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			prec += Sum.get(i,j).abs_squared();
		}
	}

	dA[site] = prec;
}


__global__ void averageGaugeQuality( double* dGff, double* dA )
{
	double gff = 0;
	double A = 0;
	for( int i = 0; i < Nx*Ny*Nz; i++ )
	{
		gff+= dGff[i];
		A += dA[i];
	}
	
	dGff[0] = gff/double(Nx*Ny*Nz)/4./3.;
	dA[0]   = A/double(Nx*Ny*Nz)/3.;
}

void generateGaugeQuality( Real** dU, int* theDevice, int deviceCount, Real**halo, Real** dHalo, double** dGff, double** dA, lat_index_t** dNnt, double &_gff, double &_A )
{
	SiteCoord<4,true> s(size);
	
	// host memory (two fields to collect the results)
	static double gff[2][32], A[2][32];
	
	for( int device=0; device<deviceCount; device++ )
	{
		gff[1][device]=0.0; 
		A[1][device]  =0.0;
	}
	for( int t=0; t<Nt; t++ )
	{
		int tDw = (t > 0)?(t-1):(s.size[0]-1);
		// do we have to copy the halo?
		if( theDevice[t] != theDevice[tDw] )
		{
			// copy the halo device -> host -> device
			cudaSetDevice(theDevice[tDw]);
			cudaMemcpy( halo[theDevice[t]], dU[tDw], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );
			cudaDeviceSynchronize();
			cudaSetDevice(theDevice[t]);
			cudaMemcpy( dHalo[theDevice[t]], halo[theDevice[t]], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
			cudaDeviceSynchronize();
			// now call kernel with dU[tDw] -> dHalo[theDevice[t]]
			generateGaugeQualityPerSite<<<Nx*Nx*Nx/32,32>>>( dU[t], dHalo[theDevice[t]], dNnt[theDevice[t]], dGff[theDevice[t]], dA[theDevice[t]] );
		}
		else
		{
			cudaSetDevice(theDevice[t]);
			generateGaugeQualityPerSite<<<Nx*Nx*Nx/32,32>>>( dU[t], dU[tDw], dNnt[theDevice[t]], dGff[theDevice[t]], dA[theDevice[t]] );
		}
		averageGaugeQuality<<<1,1>>>( dGff[theDevice[t]], dA[theDevice[t]] );
		cudaMemcpy( &gff[0][theDevice[t]], dGff[theDevice[t]], sizeof(double), cudaMemcpyDeviceToHost );
		cudaMemcpy( &A[0][theDevice[t]], dA[theDevice[t]],     sizeof(double), cudaMemcpyDeviceToHost );
		gff[1][theDevice[t]]+=gff[0][theDevice[t]];
		A[1][theDevice[t]]  +=A[0][theDevice[t]];
	}
	cudaDeviceSynchronize();
	
	for( int device=1; device<deviceCount; device++ )
	{
		gff[1][0] += gff[1][device];
		A[1][0]   += A[1][device];
	}
	_gff = gff[1][0]/(double)Nt;
	_A = A[1][0]/(double)Nt;
}

__global__ void set_hot( Real* U, int counter)
{
	typedef GpuLandauPatternParity< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;

	const lat_coord_t size[Ndim] = {1,Nx,Ny,Nz};
	SiteIndex<4,true> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;
	s.setLatticeIndex( site );
	
	PhiloxWrapper rng( site, 123, counter );

	Quaternion<Real> q;
	
	for( int mu = 0; mu < 4; mu++ )
	{
		TLinkIndex linkUp( U, s, mu );
		SU3<TLinkIndex> globUp( linkUp );

		Matrix<Complex<Real>,Nc> locMat;
		SU3<Matrix<Complex<Real>,Nc> > locU(locMat);

		locU.identity();
		
		for( int i=0; i<2; i++ )
			for( int j=i+1; j<3; j++ )
			{
				q[0] = rng.rand()*2.0-1.0;
				q[1] = rng.rand()*2.0-1.0;
				q[2] = rng.rand()*2.0-1.0;
				q[3] = rng.rand()*2.0-1.0;
				
				q.projectSU2();
				locU.rightSubgroupMult( i, j, &q );
			}
 		globUp = locU;
	}
}


int main(int argc, char* argv[])
{
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


	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	deviceCount = ( deviceCount<numbDevices ? deviceCount : numbDevices );
	
	int device;
	for( device=0; device<deviceCount; device++ ) 
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, device);
		printf("\nDevice %d: \"%s\"\n", device, deviceProp.name);
		printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);
	}
	
	
	// array to distribute the timeslices to the devices
	int theDevice[Nt];
	for( int k=0; k<deviceCount; k++ )
		for( int t = k*Nt/deviceCount; t < (k+1)*Nt/deviceCount; t++ )
			theDevice[t] = k;

// 	for( int t=0; t<Nt; t++ )
// 		cout << t << ' ' << theDevice[t] << endl;
			
	// how many timeslices does each device have?
	int numbTimeSlices[32]={0};
	for( int t=0; t<Nt; t++ )
		numbTimeSlices[theDevice[t]]++;
	
	cout << "\nTime-slice distribution:" << endl;
	cout << "Device" << '\t' << "#slices" << endl;
 	for( device=0; device<deviceCount; device++ )
 		cout << device << '\t' << numbTimeSlices[device] << endl;
    
    // min. and max. timeslice of each device
    int minT[deviceCount];
    int maxT[deviceCount];
    int midT[deviceCount];
    
    for( int t=0; t<Nt; t++ )
    {
        int tDw = (t > 0)?(t-1):(Nt-1);
        int tUp = (t == Nt-1)?0:t+1;
        
        if( theDevice[t] != theDevice[tDw] )
            minT[theDevice[t]] = t;
        else if( theDevice[t] != theDevice[tUp] )
            maxT[theDevice[t]] = t;
    }
    
    for( device=0; device<deviceCount; device++ )
    {
        midT[device] = ( minT[device] + maxT[device] ) / 2;
    }
    
    //test
    cout << "device\tminT\tmaxT\tmidT" << endl;
    for( device=0; device<deviceCount; device++ )
    {
        cout << device << '\t' << minT[device] << '\t' << maxT[device] << '\t' << midT[device] << endl;
    }
	
	// cudaStreams
	cudaStream_t streamStd[deviceCount];
	cudaStream_t streamCpy_1[deviceCount];
	cudaStream_t streamCpy_2[deviceCount];
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaStreamCreate(&streamStd[device]);
		cudaStreamCreate(&streamCpy_1[device]);
		cudaStreamCreate(&streamCpy_2[device]);
	}
	
	// cudaEvents
	cudaEvent_t eventCpyD2H[deviceCount];
	cudaEvent_t eventCpyH2D[deviceCount];
	cudaEvent_t hitFirstTimeslice[deviceCount];
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaEventCreate(&eventCpyD2H[device]);
		cudaEventCreate(&eventCpyH2D[device]);
		cudaEventCreate(&hitFirstTimeslice[device]);
	}

	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,true> s(size);
	
	// TODO maybe we should choose the filetype on compile time
// 	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,true> > lfHeaderOnly;
// 	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,true> > lfVogt;
// 	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,true> > lfPlain;


	// allocate Memory
	// host memory for configuration
// 	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for all timeslices
	Real* dU[Nt];
	for( int t=0; t<Nt; t++ )
	{
		cudaSetDevice(theDevice[t]);
		cudaMalloc( &dU[t], timesliceArraySize*sizeof(Real) );
	}
	
	// halo exchange size (only mu=0 and row 1 and 2 of SU(3))
	size_t haloSize = timesliceArraySize*sizeof(Real);
		
	// page-locked host memory for halo timeslices (one per device)
	Real* halo[deviceCount];
	for( device=0; device<deviceCount; device++ )
	{
// 		halo[device] = (Real*)malloc( haloSize );
		cudaHostAlloc( &halo[device], haloSize, 0 );
	}
	
	// device memory for halo timeslice (one per device)
	Real* dHalo[deviceCount];
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaMalloc( &dHalo[device], haloSize );
	}
	
	// device memory for collecting the parts of the gauge fixing functional and divA
	double *dGff[32];
	double *dA[32];
	
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaMalloc( &dGff[device], s.getLatticeSizeTimeslice()*sizeof(double) );
		cudaMalloc( &dA[device], s.getLatticeSizeTimeslice()*sizeof(double) );
	}

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt[32];
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaMalloc( &dNnt[device], s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t) );
	}

	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	
	// copy neighbour table to device
	for( device=0; device<deviceCount; device++ )
	{
		cudaSetDevice(device);
		cudaMemcpy( dNnt[device], nnt, s.getLatticeSizeTimeslice()*(2*(Ndim))*sizeof(lat_index_t), cudaMemcpyHostToDevice );
	}

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
// 			cudaSetDevice(theDevice[t]);
// 			cudaMemcpy( dU[t], &U[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyHostToDevice );
// 		}		
		
		
		// don't load file, fill with random numbers - then projectSU3
		cout << "\nWe fill dU with random SU(3) matrices.\n" << endl;
		int unsigned counter = 0;//time(NULL);
		
		for( int t=0; t<Nt; t++ )
		{
 			cudaSetDevice(theDevice[t]);
			set_hot<<<numBlocks*2,32>>>( dU[t], counter ); 
			counter++;
		}
		cudaDeviceSynchronize();

		// calculate and print the gauge quality
		double gff,A;
		printf( "i:\t\tgff:\t\tdA:\n");
		generateGaugeQuality( dU, theDevice, deviceCount, halo, dHalo, dGff, dA, dNnt, gff, A );
		printf( "-\t\t%1.10f\t\t%e\n", gff, A );
		
		Chronotimer kernelTimer;
		kernelTimer.reset();
		kernelTimer.start();
		
		for( int j = 0; j < orMaxIter; j++ )
		{
			for( int parity=0; parity<2; parity++ )
			{
// 				int opp_parity = parity?0:1;
				int p_offset = parity?timesliceArraySize/2:0;
				
				for( int t=0; t<Nt; t++ )
				{
					int tDw = (t > 0)?(t-1):(s.size[0]-1);
					// do we have to copy the halo?
					
					if( theDevice[t] != theDevice[tDw] )
					{
						// copy the halo: theDevice[tDw] -> host
						cudaSetDevice(theDevice[tDw]);
						cudaMemcpyAsync( halo[theDevice[t]]+p_offset, dU[tDw]+p_offset, haloSize/12, cudaMemcpyDeviceToHost, streamCpy_1[theDevice[tDw]] );
						cudaEventRecord( eventCpyD2H[theDevice[tDw]], streamCpy_1[theDevice[tDw]] );
						
						// copy the halo: host -> theDevice[t]
						cudaSetDevice(theDevice[t]);
						cudaEventSynchronize( eventCpyD2H[theDevice[tDw]] );
						cudaMemcpyAsync( dHalo[theDevice[t]]+p_offset, halo[theDevice[t]]+p_offset, haloSize/12, cudaMemcpyHostToDevice, streamCpy_2[theDevice[t]] );
						
						// now call kernel with dU[tDw] replaced by dHalo[theDevice[t]]
						orStep<<<numBlocks,threadsPerBlock,0,streamCpy_2[theDevice[t]]>>>( dU[t], dHalo[theDevice[t]], dNnt[theDevice[t]], parity, orParameter );
						cudaEventRecord( hitFirstTimeslice[theDevice[t]], streamCpy_2[theDevice[t]] );
						
						// copy halo back: theDevice[t] -> host
						cudaMemcpyAsync( halo[theDevice[t]]+p_offset, dHalo[theDevice[t]]+p_offset, haloSize/12, cudaMemcpyDeviceToHost, streamCpy_2[theDevice[t]] );
						cudaEventRecord( eventCpyD2H[theDevice[t]], streamCpy_2[theDevice[t]] );
						
						// copy halo back: host -> theDevice[tDw]
						cudaSetDevice(theDevice[tDw]);
						cudaEventSynchronize( eventCpyD2H[theDevice[t]] );
						cudaMemcpyAsync( dU[tDw]+p_offset, halo[theDevice[t]]+p_offset, haloSize/12, cudaMemcpyHostToDevice, streamCpy_2[theDevice[tDw]] );
					}
					else
					{
						cudaSetDevice(theDevice[t]);
						if( deviceCount > 1 && t == midT[theDevice[t]] )
							cudaEventSynchronize( hitFirstTimeslice[theDevice[t]] );
						orStep<<<numBlocks,threadsPerBlock,0,streamStd[theDevice[t]]>>>( dU[t], dU[tDw], dNnt[theDevice[t]], parity, orParameter );
					}
				}
				cudaDeviceSynchronize();
			}

			// calculate and print the gauge quality
			if( j % orCheckPrec == 0 )
			{
				generateGaugeQuality( dU, theDevice, deviceCount, halo, dHalo, dGff, dA, dNnt, gff, A );
				printf( "%d\t\t%1.10f\t\t%e\n", j, gff, A );
				if( A < orPrecision ) break;
			}

			totalStepNumber++;
		}
		cudaDeviceSynchronize();
		kernelTimer.stop();
		cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		
// 		cout << "Polyakov loop: " << polBefore << " - " << calculatePolyakovLoopAverage( U ) << endl;

		// copy back all timeslices
// 		for( int t=0; t<Nt; t++ )
// 			cudaMemcpy( &U[t*timesliceArraySize], dU[t], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToHost );

	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
	cout << "total kernel time: " << totalKernelTime << " s" << endl;
	cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;



}
