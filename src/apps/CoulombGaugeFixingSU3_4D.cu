/*
 * test_gaugefixing.cpp
 *
 *  Created on: Apr 18, 2012
 *      Author: vogt
 */

#include <iostream>
#include <math.h>
#include <sstream>
#include <malloc.h>
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
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
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
int fileNumberformat;
string configFile;
bool noRandomTrafo;
FileType fileType;
ReinterpretReal reinterpretReal;

// lattice setup
const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef GpuCoulombPattern<SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
typedef StandardPattern<SiteCoord<Ndim,false>,Ndim,Nc> Standard;
typedef GpuLandauPattern< SiteCoord<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice;


typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
typedef Link<GpuTimeslice,SiteCoord<Ndim-1,true>,Ndim-1,Nc> TLink3;

void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.calculateNeighbourTable( nnt );
}


//__device__ inline Real cuFabs( Real a )
//{
//	return (a>0)?(a):(-a);
//}

//__global__ void printGaugeQuality( Real* dGff, Real* dA )
//{
//	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
//	SiteCoord<3,true> s(size);
//
//	Real gff = 0;
//	Real temp = 0;
//	for( int i = 0; i < s.getLatticeSize(); i++ )
//	{
//		gff+= dGff[i];
//		if( cuFabs(dA[i]) > temp ) temp = cuFabs(dA[i]);
//	}
//
//	printf( "gff: %E\t\tdA: %E\n", gff/Real(s.getLatticeSize())/3./3., temp );
//
//}

__global__ void projectSU3( Real* U )
{
	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteCoord<3,true> s(size);
	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink3 linkUp( U, s, mu );
		SU3<TLink3> globUp( linkUp );

		globUp.projectSU3();
	}
}

//__global__ void generateGaugeQuality( Real *U, Real *dGff, Real *dA )
//{
//	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
//	SiteCoord<3,true> s(size);
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	Matrix<complex,Nc> locMatSum;
//	SU3<Matrix<complex,Nc> > Sum(locMatSum);
//
//	Sum.zero();
//
//	// TODO check if there is a faster way to compute DELTA
//	for( int mu = 1; mu < 4; mu++ )
//	{
//		s.setLatticeIndex( site );
//
//		Matrix<complex,Nc> locMat;
//		SU3<Matrix<complex,Nc> > temp(locMat);
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
//	Matrix<complex,Nc> locMatSumHerm;
//	SU3<Matrix<complex,Nc> > SumHerm(locMatSumHerm);
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
//	Matrix<complex,Nc> locTemp;
//	SU3<Matrix<complex,Nc> > temp(locTemp);
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



__global__ void __launch_bounds__(256,4) orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
{
	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.nn = nnt;

	const bool updown = threadIdx.x / 128;
	const short mu = (threadIdx.x % 128) / 32;
	const short id = (threadIdx.x % 128) % 32;

	int site = blockIdx.x * blockDim.x/8 + id;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	s.setLatticeIndex( site );
	if( (mu!=0)&&(updown==1) )
	{
		s.setNeighbour(mu-1,0);
	}

	Matrix<complex,Nc> locMat;
	SU3<Matrix<complex,Nc> > locU(locMat);

	TLink3_2 link( ((mu==0)&&(updown==1))?(UtDw):(UtUp), s, mu );

	SU3<TLink3_2> globU( link );

	// make link local
	locU.assignWithoutThirdLine(globU);
	locU.reconstructThirdLine();

	// define the update algorithm
	OrUpdate overrelax( orParameter );
	GaugeFixingSubgroupStep<SU3<Matrix<complex,Nc> >, OrUpdate, COULOMB> subgroupStep( &locU, overrelax, id, mu, updown );

	// do the subgroup iteration
	SU3<Matrix<complex,Nc> >::perSubgroup( subgroupStep );

	// copy link back
	globU.assignWithoutThirdLine(locU);
}




Real calculatePolyakovLoopAverage( Real *U )
{
	Matrix<complex,3> tempMat;
	SU3<Matrix<complex,3> > temp( tempMat );
	Matrix<complex,3> temp2Mat;
	SU3<Matrix<complex,3> > temp2( temp2Mat );

	SiteCoord<Ndim,true> s( size );

	complex result(0,0);

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
	// read parameters (command line or given config file)
	options_desc.add_options()
		("help", "produce help message")
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
		("ending", boost::program_options::value<string>(&fileEnding)->default_value(".vogt"), "file ending to append to basename (default: .vogt)")
		("postfixlabel", boost::program_options::value<string>(&postFixLabel)->default_value("_Landau"), "label to append to basename after fixing the gauge and before storing it (default _Coulomb)")
		("basename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
		("startnumber", boost::program_options::value<int>(&fileStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
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





	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	printf("\nDevice %d: \"%s\"\n", 0, deviceProp.name);
	printf("CUDA Capability Major/Minor version number:    %d.%d\n\n", deviceProp.major, deviceProp.minor);

	Chronotimer allTimer;
	allTimer.reset();

	SiteCoord<4,true> s(size);

	// TODO maybe we should choose the filetype on compile time
	LinkFile<FileHeaderOnly, Standard, Gpu, SiteCoord<4,true> > lfHeaderOnly;
	lfHeaderOnly.reinterpret = reinterpretReal;
	LinkFile<FileVogt, Standard, Gpu, SiteCoord<4,true> > lfVogt;
	lfVogt.reinterpret = reinterpretReal;
	LinkFile<FilePlain, Standard, Gpu, SiteCoord<4,true> > lfPlain;
	lfPlain.reinterpret = reinterpretReal;


	// allocate Memory
	// host memory for configuration
	Real* U = (Real*)malloc( arraySize*sizeof(Real) );

	// device memory for timeslice t
	Real* dUtUp;
	cudaMalloc( &dUtUp, timesliceArraySize*sizeof(Real) );

	// device memory for timeslice t-1
	Real* dUtDw;
	cudaMalloc( &dUtDw, timesliceArraySize*sizeof(Real) );

	// device memory for collecting the parts of the gauge fixing functional and divA
	Real *dGff;
	cudaMalloc( &dGff, s.getLatticeSizeTimeslice()*sizeof(Real) );
	Real *dA;
	cudaMalloc( &dA, s.getLatticeSizeTimeslice()*sizeof(Real) );

	// host memory for the timeslice neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt;
	cudaMalloc( &dNnt, s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof( lat_index_t ) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	// copy neighbour table to device
	cudaMemcpy( dNnt, nnt, s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = 32*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSizeTimeslice()/2/32; // // half of the lattice sites (a parity) are updated in a kernel call

	allTimer.start();

	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );

	lat_coord_t *pointerToSize;
	cudaGetSymbolAddress( (void**)&pointerToSize, "dSize" );
	GaugeFixingStats<Ndim-1,Nc,COULOMB,AVERAGE> gaugeStats( dUtUp, &size[1], pointerToSize );

	double totalKernelTime = 0;
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


			float orParameter = 1.7;

			Chronotimer kernelTimer;
			kernelTimer.reset();
			kernelTimer.start();
			for( int i = 0; i < orMaxIter; i++ )
			{
				orStep<<<numBlocks,threadsPerBlock>>>(dUtUp, dUtDw, dNnt, 0, orParameter );
				orStep<<<numBlocks,threadsPerBlock>>>(dUtUp, dUtDw, dNnt, 1, orParameter );

				if( i % orCheckPrec == 0 )
				{
//					projectSU3<<<numBlocks*2,32>>>( dUtUp );
//					projectSU3<<<numBlocks*2,32>>>( dUtDw );

					gaugeStats.generateGaugeQuality();
					printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					if( gaugeStats.getCurrentA() < orPrecision ) break;
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
