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
//#include "../lattice/gaugefixing/overrelaxation/SrUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/SU3.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/LinkFile.hxx"
//#include "../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/filetype_typedefs.h"
#include "../lattice/gaugefixing/GlobalConstants.hxx"
#include "../lattice/gaugefixing/LandauKernelsSU3.hxx"
#include "program_options/ProgramOptions.hxx"
#include "program_options/FileIterator.hxx"

using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 3;


// boost program options setup
//boost::program_options::variables_map options_vm;
//boost::program_options::options_description options_desc("Allowed options");

// parameters from command line or config file
//int nconf;
//int devicenumb;
//long seed; // TODO check datatype
//int orMaxIter;
//int orCheckPrec;
//float orParameter;
//float orPrecision;
//int reproject;
//int saSteps;
//float saMin;
//float saMax;
//int gaugeCopies;
//string fileEnding;
//string postFixLabel;
//string fileBasename;
//int fileStartnumber;
//int fileStepsize;
//int fileNumberformat;
//string configFile;
//bool noRandomTrafo;
//FileType fileType;
//ReinterpretReal reinterpretReal;

// lattice setup
//const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,NO_SPLIT>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> Gpu;


typedef Link<Gpu,SiteCoord<Ndim,FULL_SPLIT>,Ndim,Nc> TLink;


void initNeighbourTable( lat_index_t* nnt )
{
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteIndex<4,FULL_SPLIT> s( HOST_CONSTANTS::SIZE);
	s.calculateNeighbourTable( nnt );
}


__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
	SiteCoord<4,FULL_SPLIT> s(size);
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
//
//
//__global__ void __launch_bounds__(256,4) orStep( Real* U, lat_index_t* nn, bool parity, float orParameter, int counter=0  )
//{
//	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
//	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
//
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteIndex<4,FULL_SPLIT> s(size);
//	s.nn = nn;
//
//	const bool updown = threadIdx.x / NSB4;
//	const short mu = (threadIdx.x % NSB4) / NSB;
//	const short id = (threadIdx.x % NSB4) % NSB;
//
//	int site = blockIdx.x * blockDim.x/8 + id;
//	if( parity == 1 ) site += s.getLatticeSize()/2;
//
//	s.setLatticeIndex( site );
//	if( updown==1 )
//	{
//		s.setNeighbour(mu,false);
//	}
//
////	if(id == 0) printf("bin in or\n");
//
//	Matrix<Complex<Real>,Nc> locMat;
//	SU3<Matrix<Complex<Real>,Nc> > locU(locMat);
//
//	TLinkIndex link( U, s, mu );
//
//	SU3<TLinkIndex> globU( link );
//
//	// make link local
//	locU.assignWithoutThirdLine(globU);
//	locU.reconstructThirdLine();
//
//	// define the update algorithm
//	OrUpdate overrelax( orParameter );
//	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, OrUpdate, LANDAU> subgroupStep( &locU, overrelax, id, mu, updown );
////	MicroUpdate micro;
////	GaugeFixingSubgroupStep<SU3<Matrix<Complex<Real>,Nc> >, MicroUpdate, LANDAU> subgroupStep( &locU, micro, id, mu, updown );
//
//	// do the subgroup iteration
//	SU3<Matrix<Complex<Real>,Nc> >::perSubgroup( subgroupStep );
//
//	// copy link back
////	globU=locU; //TODO with or without 3rd line?
//	globU.assignWithoutThirdLine(locU);
//
//	// project back
////	globU.projectSU3withoutThirdRow();
//}
//
//// TODO Hack, because cuda5.0 has a bug and the native *= operator does not work
//__device__ Matrix<Complex<Real>,3>& mult( Matrix<Complex<Real>,3>& c, Matrix<Complex<Real>,3>& a, Matrix<Complex<Real>,3>& b )
//{
//	for(int i = 0; i < 3; i++ )
//	{
//		for( int j = 0; j < 3; j++ )
//		{
//			c.mat[i*3+j] = 0;
//			for( int k = 0; k < 3; k++ )
//			{
//				c.mat[i*3+j] += a.mat[i*3+k] * b.mat[k*3+j];
//			}
//		}
//	}
//	return c;
//}
//
//__global__ void calculatePlaquette( Real *U, lat_index_t* nn, double *dPlaquette )
//{
//	typedef GpuLandauPattern< SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> GpuIndex;
//	typedef Link<GpuIndex,SiteIndex<Ndim,FULL_SPLIT>,Ndim,Nc> TLinkIndex;
//
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteIndex<4,FULL_SPLIT> s(size);
//	s.nn = nn;
//
//	Matrix<Complex<Real>,Nc> matP;
////	SU3<Matrix<Complex<Real>,Nc> > P(matP);
//
//	Matrix<Complex<Real>,Nc> matTemp;
//	SU3<Matrix<Complex<Real>,Nc> > temp(matTemp);
//
//	Matrix<Complex<Real>,Nc> matTemp2;
//
//
//	double localPlaquette = 0;
//
//
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		for( int nu = mu+1; nu < 4; nu++)
//		{
//			for( int i = 0; i < 3; i++ )
//				for( int j =0; j < 3; j++ )
//				{
//					if( i == j ) matP.set( i,j,Complex<Real>(1.0,.0) );
//					else matP.set( i,j,Complex<Real>(.0,.0) );
//				}
//
//			{
//				s.setLatticeIndex( site );
//
//				TLinkIndex link( U, s, mu );
//				SU3<TLinkIndex> globU( link );
//
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//
//				mult(matTemp2, matP, temp.mat );
//				matP = matTemp2;
////				matTemp2 = matP * matTemp;
//			}
//
//			{
//				s.setNeighbour( mu, true );
//
//				TLinkIndex link( U, s, nu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//			{
//				s.setLatticeIndex( site );
//				s.setNeighbour(nu, true );
//
//				TLinkIndex link( U, s, mu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//				temp.hermitian();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//			{
//				s.setLatticeIndex( site );
//
//				TLinkIndex link( U, s, nu );
//				SU3<TLinkIndex> globU( link );
//				temp.assignWithoutThirdLine( globU );
//				temp.reconstructThirdLine();
//				temp.hermitian();
//
//				mult(matTemp2,matP, temp.mat );
//				matP = matTemp2;
//			}
//
//
//
//
//
//			localPlaquette += matP.trace().x;
//		}
//	}
//
//	dPlaquette[site] = localPlaquette/6./3.;
////	dPlaquette[site] = 1;
////	dPlaquette[site] = 0;
//}
//
//__global__ void printPlaquette( double* dPlaquette )
//{
//	const lat_coord_t size[Ndim] = {Nx,Ny,Nz,Nt};
//	SiteCoord<4,FULL_SPLIT> s(size);
//
//	double plaquette = 0;
//	for( int i = 0; i < s.getLatticeSize(); i++ )
//	{
//		plaquette += dPlaquette[i];
//	}
//
//	printf( "\t%E\n", plaquette/double(s.getLatticeSize()) );
//
//}

//Real calculatePolyakovLoopAverage( Real *U )
//{
//	Matrix<Complex<Real>,3> tempMat;
//	SU3<Matrix<Complex<Real>,3> > temp( tempMat );
//	Matrix<Complex<Real>,3> temp2Mat;
//	SU3<Matrix<Complex<Real>,3> > temp2( temp2Mat );
//
//	SiteCoord<Ndim,FULL_SPLIT> s( HOST_CONSTANTS::SIZE );
//
//	Complex<Real> result(0,0);
//
//	for( s[1] = 0; s[1] < s.size[1]; s[1]++ )
//	{
//		for( s[2] = 0; s[2] < s.size[2]; s[2]++ )
//		{
//			for( s[3] = 0; s[3] < s.size[3]; s[3]++ )
//			{
//				temp.identity();
//				temp2.zero();
//
//				for( s[0] = 0; s[0] < s.size[0]; s[0]++ )
//				{
//
//					TLink link( U, s, 0 );
//					SU3<TLink> globU( link );
//
//					temp2 = temp2 + temp*globU;
//
//					temp = temp2;
//					temp2.zero();
//				}
//				result += temp.trace();
//			}
//		}
//	}
//
//	return sqrt(result.x*result.x+result.y*result.y) / (Real)(s.getLatticeSizeTimeslice()*Nc);
//}






int main(int argc, char* argv[])
{
	LandauKernelsSU3::initCacheConfig();

	// read configuration from file or command line
	ProgramOptions options;
	int returncode = options.init( argc, argv );
	if( returncode != 0 ) return returncode;


	// read parameters (command line or given config file)
//	options_desc.add_options()
//		("help", "produce help message")
//		("devicenumb", boost::program_options::value<int>(&devicenumb)->default_value(0), "which CUDA device to use")
//		("nconf,m", boost::program_options::value<int>(&nconf)->default_value(1), "how many files to gaugefix")
//		("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(1000), "Max. number of OR iterations")
//		("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")
//		("sasteps", boost::program_options::value<int>(&saSteps)->default_value(1000), "number of SA steps")
//		("samin", boost::program_options::value<float>(&saMin)->default_value(.01), "min. SA temperature")
//		("samax", boost::program_options::value<float>(&saMax)->default_value(.4), "max. SA temperature")
//		("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "OR parameter")
//		("orprecision", boost::program_options::value<float>(&orPrecision)->default_value(1E-7), "OR precision (dmuAmu)")
//		("orcheckprecision", boost::program_options::value<int>(&orCheckPrec)->default_value(100), "how often to check the gauge precision")
//		("reproject", boost::program_options::value<int>(&reproject)->default_value(100), "reproject every arg-th step")
//		("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "Number of gauge copies")
//		("ending", boost::program_options::value<string>(&fileEnding)->default_value(".vogt"), "file ending to append to basename")
//		("postfixlabel", boost::program_options::value<string>(&postFixLabel)->default_value("_Landau"), "label to append to basename after fixing the gauge and before storing it")
//		("basename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
//		("startnumber", boost::program_options::value<int>(&fileStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
//		("stepsize", boost::program_options::value<int>(&fileStepsize)->default_value(1), "file numbering startnumber, startnumber+stepsize,...")
//		("numberformat", boost::program_options::value<int>(&fileNumberformat)->default_value(1), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
//		("filetype", boost::program_options::value<FileType>(&fileType), "type of configuration (PLAIN, HEADERONLY, VOGT)")
//		("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")
//		("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpretReal)->default_value(STANDARD), "reinterpret Real datatype (STANDARD = do nothing, FLOAT = convert input as float and cast to Real, DOUBLE = ...)")
//		("norandomtrafo", boost::program_options::value<bool>(&noRandomTrafo)->default_value(false), "no random gauge trafo" )
//		;
//
//	boost::program_options::positional_options_description options_p;
//	options_p.add("config-file", -1);
//
//	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
//			options(options_desc).positional(options_p).run(), options_vm);
//	boost::program_options::notify(options_vm);
//
//	ifstream cfg( configFile.c_str() );
//	boost::program_options::store(boost::program_options::parse_config_file( cfg, options_desc), options_vm);
//	boost::program_options::notify(options_vm);
//
//	if (options_vm.count("help")) {
//		cout << "Usage: " << argv[0] << " [options] [config-file]" << endl;
//		cout << options_desc << "\n";
//		return 1;
//	}



	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, options.getDeviceNumber() );
	cudaSetDevice( options.getDeviceNumber() );

	printf("\nDevice %d: \"%s\"\n", options.getDeviceNumber(), deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);

	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,FULL_SPLIT> s(HOST_CONSTANTS::SIZE);


	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfHeaderOnly;
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfVogt;
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,FULL_SPLIT> > lfPlain;


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


	double* dPlaquette;
	cudaMalloc( &dPlaquette, s.getLatticeSize()*sizeof(double) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nn );
	// copy neighbour table to device
	cudaMemcpy( dNn, nn, s.getLatticeSize()*(2*(Ndim))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = NSB*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSize()/2/NSB; // // half of the lattice sites (a parity) are updated in a kernel call

	allTimer.start();

	GaugeFixingStats<Ndim,Nc,LandauKernelsSU3,AVERAGE> gaugeStats( dU, HOST_CONSTANTS::SIZE );


	double totalKernelTime = 0;

	long totalStepNumber = 0;

	FileIterator fi( options );
	for( fi.reset(); fi.hasNext(); fi.next() )
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

		// copying configuration ...
		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

		// calculate and print the gauge quality
		printf( "i:\t\tgff:\t\tdA:\n");
		gaugeStats.generateGaugeQuality();
		printf( "   \t\t%1.10f\t\t%e\n", gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

		Chronotimer kernelTimer;
		kernelTimer.reset();
		kernelTimer.start();



		float temperature = options.getSaMax();
		float tempStep = (options.getSaMax()-options.getSaMin())/(float)options.getSaSteps();

		for( int i = 0; i < options.getSaSteps(); i++ )
		{


			LandauKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 0, temperature, PhiloxWrapper::getNextCounter() );
			LandauKernelsSU3::saStep(numBlocks,threadsPerBlock,dU, dNn, 1, temperature, PhiloxWrapper::getNextCounter() );

			for( int mic = 0; mic < options.getSaMicroupdates(); mic++ )
			{
				LandauKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 0 );
				LandauKernelsSU3::microStep(numBlocks,threadsPerBlock,dU, dNn, 1 );
			}


			if( i % options.getOrCheckPrecision() == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );


				gaugeStats.generateGaugeQuality();
//				CudaError::getLastError( "generateGaugeQuality error" );
				printf( "%d\t%f\t\t%1.10f\t\t%e\n", 0, temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );
			}
			temperature -= tempStep;
		}


		for( int i = 0; i < options.getOrMaxIter(); i++ )
		{

			LandauKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 0, options.getOrParameter() );
			LandauKernelsSU3::orStep(numBlocks,threadsPerBlock,dU, dNn, 1, options.getOrParameter() );

			if( i % options.getOrCheckPrecision() == 0 )
			{
				projectSU3<<<numBlocks*2,32>>>( dU );
				gaugeStats.generateGaugeQuality();
				printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

				if( gaugeStats.getCurrentA() < options.getOrPrecision() ) break;
			}

			totalStepNumber++;
		}

		cudaThreadSynchronize();
		kernelTimer.stop();
		cout << "kernel time: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();

		// copy back
		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );
		
		//saving file
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

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
	cout << "total kernel time: " << totalKernelTime << " s" << endl;

	cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;

}
