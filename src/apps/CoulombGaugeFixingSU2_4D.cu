/*
 *
 *  Created on: Apr 18, 2012
 *      Author: vogt
 */

#include <iostream>
#include <math.h>
#include <sstream>
#include <malloc.h>
#include "../lattice/gaugefixing/GaugeFixingStepSU2.hxx"
#include "../lattice/gaugefixing/GaugeFixingStats.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrUpdate.hxx"
#include "../lattice/gaugefixing/overrelaxation/MicroUpdate.hxx"
#include "../lattice/gaugefixing/simulated_annealing/SaUpdate.hxx"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"
#include "../lattice/Matrix.hxx"
#include "../lattice/SU2.hxx"
#include "../lattice/LinkFile.hxx"
#include "../lattice/gaugefixing/overrelaxation/OrSubgroupStep.hxx"
#include "../util/timer/Chronotimer.h"
#include "../lattice/filetypes/FileHeaderOnly.hxx"
#include "../lattice/filetypes/FilePlain.hxx"
#include "../lattice/filetypes/FileVogt.hxx"
#include "../lattice/filetypes/filetype_typedefs.h"
#include "../util/rng/PhiloxWrapper.hxx"
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include "../lattice/gaugefixing/GlobalConstants.hxx"


using namespace std;

const lat_dim_t Ndim = 4;
const short Nc = 2;

//#ifdef _X_
//const lat_coord_t Nx = _X_;
//#else
//#error "Define X (the lattice size in x-direction)"
//#endif
//#ifdef _Y_
//const lat_coord_t Ny = _Y_;
//#else
//const lat_coord_t Ny = _X_;
//bool warnY = true; // TODO print the warning
//#endif
//#ifdef _Z_
//const lat_coord_t Nz = _Z_;
//#else
//const lat_coord_t Nz = _X_;
//bool warnZ = true;
//#endif
//#ifdef _T_
//const lat_coord_t Nt = _T_;
//#else
//#error "Define T (the lattice size in t-direction)"
//#endif


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
int reproject;
int saSteps;
float saMin;
float saMax;
int microupdates;
int gaugeCopies;
string fileEnding;
string fileBasename;
int fileStartnumber;
int fileNumberformat;
string configFile;
bool noRandomTrafo;
FileType fileType;
ReinterpretReal reinterpretReal;

// lattice setup
//const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//__constant__ lat_coord_t dSize[Ndim] = {Nt,Nx,Ny,Nz};
const int arraySize = Nt*Nx*Ny*Nz*Ndim*Nc*Nc*2;
const int timesliceArraySize = Nx*Ny*Nz*Ndim*Nc*Nc*2;

typedef StandardPattern<SiteCoord<Ndim,false>,Ndim,Nc> Standard;
typedef GpuCoulombPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;


typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;


__device__ inline Real cuFabs( Real a )
{
	return (a>0)?(a):(-a);
}

void initNeighbourTable( lat_index_t* nnt )
{
	const lat_coord_t size[Ndim] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.calculateNeighbourTable( nnt );
}


__device__ inline Quaternion<Real> sumOverStaples( Real* UtDw, Real* Ut, Real* UtUp, SiteIndex<3,true> *site, lat_dim_t mu )
{
	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	Quaternion<Real> A;
	A.zero();
	Quaternion<Real> temp;
	lat_index_t saveSite = site->getLatticeIndex();
	for( int nu = 0; nu < 4; nu++ )
	{
		site->setLatticeIndex( saveSite );
		if( nu != mu )
		{
			// The following is rather complex since we have to carefully set in which timeslice we are...
			// For each link we have to set three things: the direction (mu/nu), the timeslice (UtDw, Ut, UtUp), the site (via setNeighbour)
			// and choose if we want the hermitian conjugate
			if( mu != 0 ) site->setNeighbour( mu-1, true );
			TLink3_2 link( (mu==0)?UtUp:Ut, *site, nu );
			SU2<TLink3_2> glob( link );

			temp = glob.getQuaternion();


			glob.getMat().getSite().setLatticeIndex( saveSite );
			glob.getMat().setMu( mu );
			if( nu != 0 )
			{
				glob.getMat().getSite().setNeighbour( nu-1, true );
				glob.getMat().setPointer( Ut );
			}
			else
			{
				glob.getMat().setPointer( UtUp );
			}
			temp *= glob.getQuaternionHermitian();


			glob.getMat().getSite().setLatticeIndex( saveSite );
			glob.getMat().setMu( nu );
			glob.getMat().setPointer( Ut );
			temp *= glob.getQuaternionHermitian();


			A += temp;


			glob.getMat().getSite().setLatticeIndex( saveSite );
			glob.getMat().setMu( nu );
			if( mu != 0 )
			{
				glob.getMat().getSite().setNeighbour( mu-1, true );
				if( nu != 0 )
				{
					glob.getMat().getSite().setNeighbour( nu-1, false );
					glob.getMat().setPointer( Ut );
				}
				else
				{
					glob.getMat().setPointer( UtDw );
				}
			}
			else
			{
				glob.getMat().setPointer( UtUp );
				glob.getMat().getSite().setNeighbour( nu-1, false ); // if mu == 0 then nu != 0 so nothing to check here!
			}

			temp = glob.getQuaternionHermitian();


			glob.getMat().getSite().setLatticeIndex( saveSite );
			glob.getMat().setMu( mu );
			if( nu != 0 )
			{
				glob.getMat().getSite().setNeighbour( nu-1, false );
				glob.getMat().setPointer( Ut );
			}
			else
			{
				glob.getMat().setPointer( UtDw );
			}
			temp *= glob.getQuaternionHermitian();


			glob.getMat().getSite().setLatticeIndex( saveSite );
			glob.getMat().setMu( nu );
			if( nu != 0 )
			{
				glob.getMat().getSite().setNeighbour( nu-1, false );
				glob.getMat().setPointer( Ut );
			}
			else
			{
				glob.getMat().setPointer( UtDw );
			}
			temp *= glob.getQuaternion();

			A += temp;
//			if( )
		}
	}
	return A;

}

__global__ void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, 5168, counter ); // TODO take care that we do not reuse random numbers!!!

	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.nn = nnt;

	int site = blockIdx.x * blockDim.x + threadIdx.x;
	if( parity == 1 ) site += s.getLatticeSize()/2;

	// calculate sum over staples
	for( lat_dim_t mu = 0; mu < 4; mu ++ )
	{
		Quaternion<Real> A = sumOverStaples( UtDw, Ut, UtUp, &s, mu );




		Real e0,e1,e2,e3, dk, p0;
		Real r1,r2,r3,r4;
		Real a0,a1,a2,a3;
		Real delta, phi, sin_alpha, sin_theta, cos_theta;
		e0=A[0];
		e1=A[1]; //
		e2=A[2]; //
		e3=A[3]; //
		dk=rsqrt(e0*e0+e1*e1+e2*e2+e3*e3);
		p0=(dk*1./beta);  // entspricht abeta


	//	if(beta_eff<=0.0)
	//	  cout << "beta is zero";

		do
		{
		  do; while ((r1 = rng.rand()) < 0.0001);
		  r1 = -log(r1)*p0;
		  do; while ((r2 = rng.rand()) < 0.0001);
		  r2 = -log(r2)*p0;
		  r3 = cospi(2.0*rng.rand());
		  r3 = r3*r3;
		  delta = r2+r1*r3;
		  r4=rng.rand();
		} while (r4*r4 > (1.0-0.5*delta));
		a0=1.0-delta;
		cos_theta=2.0*rng.rand()-1.0;
		sin_theta=sqrt(1.0-cos_theta*cos_theta);
		sin_alpha=sqrt(1-a0*a0);
		phi=2.0*rng.rand();
		a1=sin_alpha*sin_theta*cospi(phi);
		a2=sin_alpha*sin_theta*sinpi(phi);
		a3=sin_alpha*cos_theta;


	//	dk = 1/dk;

		e0 *= dk;
		e1 *= dk;
		e2 *= dk;
		e3 *= dk;

		Quaternion<Real> X;

		X[0] = a0*e0+a3*e3+a2*e2+e1*a1;
		X[3] = e0*a3-e3*a0+a1*e2-a2*e1;
		X[2] = a3*e1-a0*e2+a2*e0-a1*e3;
		X[1] = a2*e3+a1*e0-a3*e2-e1*a0;

		A.hermitian();
		X *= A;


		s.setLatticeIndex( site );
		TLink3_2 link( Ut, s, mu );
		SU2<TLink3_2> glob( link );

		glob.setQuaternion( X );
	}


}

__global__ void restoreSecondLine( Real* U, lat_index_t* nnt )
{
	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.nn = nnt;

	int site = blockIdx.x * blockDim.x + threadIdx.x;

	s.setLatticeIndex( site );

	for( int mu = 0; mu < 4; mu++ )
	{
		TLink3_2 link( U, s, mu );
		SU2<TLink3_2> glob( link );
		glob.projectSU2();
	}
}


__global__ void randomTrafo( Real* UtUp, Real* UtDw,lat_index_t* nnt, bool parity, int counter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, 321, counter );

	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.nn = nnt;


	int site = blockIdx.x * blockDim.x + threadIdx.x;
	if( parity == 1 ) site += s.getLatticeSize()/2;


	Quaternion<Real> randTrafo;

	Real alpha, phi, cos_theta, sin_theta, sin_alpha;
	alpha = rng.rand();
	phi = 2.0 * rng.rand();
	cos_theta = 2.0 * rng.rand() - 1.0;
	sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	sin_alpha = sinpi(alpha);
	randTrafo[0] = cospi(alpha);
	randTrafo[1] = sin_alpha * sin_theta * cospi(phi);
	randTrafo[2] = sin_alpha * sin_theta * sinpi(phi);
	randTrafo[3] = sin_alpha * cos_theta;

//		randTrafo[0] = rng.rand();
//		randTrafo[1] = rng.rand();
//		randTrafo[2] = rng.rand();
//		randTrafo[3] = rng.rand();

//		randTrafo.projectSU2();
	Quaternion<Real> randTrafoHermitian;
	randTrafoHermitian = randTrafo;
	randTrafoHermitian.hermitian();

	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		TLink3_2 linkUp( UtUp, s, mu );
		SU2<TLink3_2> globUp( linkUp );

		Quaternion<Real> locUp;

		Complex<Real> a = globUp.get(0,0);
		locUp[0] = a.x;
		locUp[3] = a.y;
		a = globUp.get(0,1);
		locUp[2] = a.x;
		locUp[1] = a.y;


		if( mu != 0 ) s.setNeighbour( mu-1, false );
		TLink3_2 linkDw( (mu != 0)?UtUp:UtDw, s, mu );
		SU2<TLink3_2> globDw( linkDw );

		Quaternion<Real> locDw;

		a = globDw.get(0,0);
		locDw[0] = a.x;
		locDw[3] = a.y;
		a = globDw.get(0,1);
		locDw[2] = a.x;
		locDw[1] = a.y;




		locUp.leftMult( randTrafo );
		locDw.rightMult( randTrafoHermitian );

		a.x = locUp[0];
		a.y = locUp[3];
		globUp.set(0,0,a);
		a.x = locUp[2];
		a.y = locUp[1];
		globUp.set(0,1,a);

		a.x = locDw[0];
		a.y = locDw[3];
		globDw.set(0,0,a);
		a.x = locDw[2];
		a.y = locDw[1];
		globDw.set(0,1,a);
	}
}

/**
 * TODO does not work...
 */
__global__ void orStepSingleThread( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
{
	typedef GpuLandauPattern< SiteIndex<Ndim-1,true>,Ndim-1,Nc> GpuTimeslice_2;
	typedef Link<GpuTimeslice_2,SiteIndex<Ndim-1,true>,Ndim-1,Nc> TLink3_2;

	const lat_coord_t size[Ndim-1] = {Nx,Ny,Nz};
	SiteIndex<3,true> s(size);
	s.nn = nnt;


	int site = blockIdx.x * blockDim.x + threadIdx.x;
	if( parity == 1 ) site += s.getLatticeSize()/2;


//	if( (mu!=0)&&(updown==1) )
//	{
//		s.setNeighbour(mu-1,0);
//	}


	// make 8 quaternions local
	Quaternion<Real> locQuat[8];

	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		int index = mu*2;

		TLink3_2 link( UtUp, s, mu );
		SU2<TLink3_2> globU( link );

		Complex<Real> a = globU.get(0,0);
		locQuat[index][0] = a.x;
		locQuat[index][3] = a.y;
		a = globU.get(0,1);
		locQuat[index][2] = a.x;
		locQuat[index][1] = a.y;
	}

	for( int mu = 1; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		s.setNeighbour(mu-1,0);
		int index = mu*2+1;

		TLink3_2 link( UtUp, s, mu );
		SU2<TLink3_2> globU( link );

		Complex<Real> a = globU.get(0,0);
		locQuat[index][0] = a.x;
		locQuat[index][3] = a.y;
		a = globU.get(0,1);
		locQuat[index][2] = a.x;
		locQuat[index][1] = a.y;
	}

	{
		s.setLatticeIndex( site );
		int index = 1;
		TLink3_2 link( UtDw, s, 0 );
		SU2<TLink3_2> globU( link );

		Complex<Real> a = globU.get(0,0);
		locQuat[index][0] = a.x;
		locQuat[index][3] = a.y;
		a = globU.get(0,1);
		locQuat[index][2] = a.x;
		locQuat[index][1] = a.y;
	}


	Quaternion<Real> A;
	A[0] = 0;
	A[1] = 0;
	A[2] = 0;
	A[3] = 0;

	for( int mu = 2; mu < 8; mu+=2 )
	{
		A[0] += locQuat[mu][0];
		A[1] -= locQuat[mu][1];
		A[2] -= locQuat[mu][2];
		A[3] -= locQuat[mu][3];
	}
	for( int mu = 3; mu < 8; mu+=2 )
	{
		A[0] += locQuat[mu][0];
		A[1] += locQuat[mu][1];
		A[2] += locQuat[mu][2];
		A[3] += locQuat[mu][3];
	}



//	if( id == 0 && mu == 0 && updown == 0 ) printf( "%f\n", globU.get(0,0).x );

	// define the update algorithm
	Real ai_sq = A[1]*A[1]+A[2]*A[2]+A[3]*A[3];
	Real a0_sq = A[0]*A[0];

	Real b=(orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
	Real c=rsqrt(a0_sq+b*b*ai_sq);

	A[0]*=c;
	A[1]*=b*c;
	A[2]*=b*c;
	A[3]*=b*c;


	for( int mu = 0; mu < 8; mu+=2 )
	{
		locQuat[mu].leftMult( A );
	}
	for( int mu = 1; mu < 8; mu+=2 )
	{
		A.hermitian();
		locQuat[mu].rightMult( A );
	}

//	globU.assignQuaternion(locU);

	// reproject
	for( int mu = 0; mu < 8; mu++ )
		locQuat[mu].projectSU2();



	for( int mu = 0; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		int index = mu*2;

		TLink3_2 link( UtUp, s, mu );
		SU2<TLink3_2> globU( link );

		Complex<Real> a;
		a.x = locQuat[index][0];
		a.y = locQuat[index][3];
		globU.set(0,0,a);

		a.x = locQuat[index][2];
		a.y = locQuat[index][1];
		globU.set(0,1,a);
	}

	for( int mu = 1; mu < 4; mu++ )
	{
		s.setLatticeIndex( site );
		s.setNeighbour(mu-1,0);
		int index = mu*2+1;

		TLink3_2 link( UtUp, s, mu );
		SU2<TLink3_2> globU( link );

		Complex<Real> a;
		a.x = locQuat[index][0];
		a.y = locQuat[index][3];
		globU.set(0,0,a);

		a.x = locQuat[index][2];
		a.y = locQuat[index][1];
		globU.set(0,1,a);
	}

	{
		s.setLatticeIndex( site );
		int index = 1;
		TLink3_2 link( UtDw, s, 0 );
		SU2<TLink3_2> globU( link );

		Complex<Real> a;
		a.x = locQuat[index][0];
		a.y = locQuat[index][3];
		globU.set(0,0,a);

		a.x = locQuat[index][2];
		a.y = locQuat[index][1];
		globU.set(0,1,a);
	}


}


__global__ void orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter )
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



//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtUp adress: %p \n", UtUp);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtDw adress: %p \n", UtDw);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("nntadress: %p \n", nnt);

//	Matrix<complex,Nc> locMat;
//	SU2<Matrix<complex,Nc> > locU(locMat);

	Quaternion<Real> locQuat;
//	SU2<Quaternion<Real> > locU(locQuat);

	TLink3_2 link( ((mu==0)&&(updown==1))?(UtDw):(UtUp), s, mu );
	SU2<TLink3_2> globU( link );

	// make link local
//	locU.assignQuaternion(globU);

	Complex<Real> a = globU.get(0,0);
	locQuat[0] = a.x;
	locQuat[3] = a.y;
	a = globU.get(0,1);
	locQuat[2] = a.x;
	locQuat[1] = a.y;


//	if( id == 0 && mu == 0 && updown == 0 ) printf( "%f\n", globU.get(0,0).x );

	// define the update algorithm
	OrUpdate overrelax( orParameter );
	GaugeFixingStepSU2<OrUpdate, COULOMB> TheUpdate( &locQuat, overrelax, id, mu, updown );

	TheUpdate.apply();



//	globU.assignQuaternion(locU);

	// reproject
	locQuat.projectSU2();

	// copy link back
	a.x = locQuat[0];
	a.y = locQuat[3];
	globU.set(0,0,a);
	a.x = locQuat[2];
	a.y = locQuat[1];
	globU.set(0,1,a);
}



__global__ void microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity )
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



//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtUp adress: %p \n", UtUp);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtDw adress: %p \n", UtDw);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("nntadress: %p \n", nnt);

//	Matrix<complex,Nc> locMat;
//	SU2<Matrix<complex,Nc> > locU(locMat);

	Quaternion<Real> locQuat;
//	SU2<Quaternion<Real> > locU(locQuat);

	TLink3_2 link( ((mu==0)&&(updown==1))?(UtDw):(UtUp), s, mu );
	SU2<TLink3_2> globU( link );

	// make link local
//	locU.assignQuaternion(globU);

	Complex<Real> a = globU.get(0,0);
	locQuat[0] = a.x;
	locQuat[3] = a.y;
	a = globU.get(0,1);
	locQuat[2] = a.x;
	locQuat[1] = a.y;


//	if( id == 0 && mu == 0 && updown == 0 ) printf( "%f\n", globU.get(0,0).x );

	// define the update algorithm
	MicroUpdate microupdate;
	GaugeFixingStepSU2<MicroUpdate, COULOMB> TheUpdate( &locQuat, microupdate, id, mu, updown );

	TheUpdate.apply();



//	globU.assignQuaternion(locU);

	// reproject
	locQuat.projectSU2();

	// copy link back
	a.x = locQuat[0];
	a.y = locQuat[3];
	globU.set(0,0,a);
	a.x = locQuat[2];
	a.y = locQuat[1];
	globU.set(0,1,a);
}





__global__ void saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int counter )
{
	PhiloxWrapper rng( blockIdx.x * blockDim.x + threadIdx.x, 123, counter );

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



//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtUp adress: %p \n", UtUp);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("UtDw adress: %p \n", UtDw);
//	if(blockIdx.x * blockDim.x + threadIdx.x == 0  ) printf("nntadress: %p \n", nnt);

//	Matrix<complex,Nc> locMat;
//	SU2<Matrix<complex,Nc> > locU(locMat);

	Quaternion<Real> locQuat;
//	SU2<Quaternion<Real> > locU(locQuat);

	TLink3_2 link( ((mu==0)&&(updown==1))?(UtDw):(UtUp), s, mu );
	SU2<TLink3_2> globU( link );

	// make link local
//	locU.assignQuaternion(globU);

	Complex<Real> a = globU.get(0,0);
	locQuat[0] = a.x;
	locQuat[3] = a.y;
	a = globU.get(0,1);
	locQuat[2] = a.x;
	locQuat[1] = a.y;


//	if( id == 0 && mu == 0 && updown == 0 ) printf( "%f\n", globU.get(0,0).x );

	// define the update algorithm
	SaUpdate sa( temperature, &rng );
	GaugeFixingStepSU2<SaUpdate, COULOMB> TheUpdate( &locQuat, sa, id, mu, updown );

	TheUpdate.apply();



//	globU.assignQuaternion(locU);

	// reproject
	locQuat.projectSU2();

	// copy link back
	a.x = locQuat[0];
	a.y = locQuat[3];
	globU.set(0,0,a);
	a.x = locQuat[2];
	a.y = locQuat[1];
	globU.set(0,1,a);
}













//__global__ void calculatePlaquette( Real *U, lat_index_t* nn, double *dPlaquette )
//{
//	typedef GpuLandauPattern< SiteIndex<Ndim,true>,Ndim,Nc> GpuIndex;
//	typedef Link<GpuIndex,SiteIndex<Ndim,true>,Ndim,Nc> TLinkIndex;
//
//	int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//	const lat_coord_t size[Ndim] = {Nt,Nx,Ny,Nz};
//	SiteIndex<4,true> s(size);
//	s.nn = nn;
//
//	Matrix<complex,Nc> matP;
//	SU3<Matrix<complex,Nc> > P(matP);
//
//	Matrix<complex,Nc> matTemp;
//	SU3<Matrix<complex,Nc> > temp(matTemp);
//
//
//	double localPlaquette = 0;
//
//
//	for( int mu = 0; mu < 4; mu++ )
//	{
//		for( int nu = mu+1; nu < 4; nu++)
//		{
//			P.identity();
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
//				P *= temp;
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
//				P *= temp;
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
//				P *= temp;
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
//				P *= temp;
//			}
//
//
//
//
//
//			localPlaquette += P.trace().x;
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
//	SiteCoord<4,true> s(size);
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
//	Matrix<complex,3> tempMat;
//	SU3<Matrix<complex,3> > temp( tempMat );
//	Matrix<complex,3> temp2Mat;
//	SU3<Matrix<complex,3> > temp2( temp2Mat );
//
//	SiteCoord<Ndim,true> s( size );
//
//	complex result(0,0);
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
	bool doMc;


	// read parameters (command line or given config file)
	options_desc.add_options()
		("help", "produce help message")
		("nconf,m", boost::program_options::value<int>(&nconf)->default_value(1), "how many files to gaugefix")
		("ormaxiter", boost::program_options::value<int>(&orMaxIter)->default_value(1000), "Max. number of OR iterations")
		("seed", boost::program_options::value<long>(&seed)->default_value(1), "RNG seed")
		("sasteps", boost::program_options::value<int>(&saSteps)->default_value(1000), "number of SA steps")
		("samin", boost::program_options::value<float>(&saMin)->default_value(.01), "min. SA temperature")
		("samax", boost::program_options::value<float>(&saMax)->default_value(.4), "max. SA temperature")
		("microupdates", boost::program_options::value<int>(&microupdates)->default_value(3), "number of microcanoncial updates at each SA temperature")
		("orparameter", boost::program_options::value<float>(&orParameter)->default_value(1.7), "OR parameter")
		("orprecision", boost::program_options::value<float>(&orPrecision)->default_value(1E-7), "OR precision (dmuAmu)")
		("orcheckprecision", boost::program_options::value<int>(&orCheckPrec)->default_value(100), "how often to check the gauge precision")
		("reproject", boost::program_options::value<int>(&reproject)->default_value(100), "reproject every arg-th step")
		("gaugecopies", boost::program_options::value<int>(&gaugeCopies)->default_value(1), "Number of gauge copies")
		("ending", boost::program_options::value<string>(&fileEnding)->default_value(".vogt"), "file ending to append to basename (default: .vogt)")
		("basename", boost::program_options::value<string>(&fileBasename), "file basename (part before numbering starts)")
		("startnumber", boost::program_options::value<int>(&fileStartnumber)->default_value(0), "file index number to start from (startnumber, ..., startnumber+nconf-1")
		("numberformat", boost::program_options::value<int>(&fileNumberformat)->default_value(1), "number format for file index: 1 = (0,1,2,...,10,11), 2 = (00,01,...), 3 = (000,001,...),...")
		("filetype", boost::program_options::value<FileType>(&fileType), "type of configuration (PLAIN, HEADERONLY, VOGT)")
		("config-file", boost::program_options::value<string>(&configFile), "config file (command line arguments overwrite config file settings)")
		("reinterpret", boost::program_options::value<ReinterpretReal>(&reinterpretReal)->default_value(STANDARD), "reinterpret Real datatype (STANDARD = do nothing, FLOAT = read input as float and cast to Real, DOUBLE = ...)")
		("domc", boost::program_options::value<bool>(&doMc)->default_value(false), "do the MC tests...")


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

	SiteCoord<4,true> s(HOST_CONSTANTS::SIZE);


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

	// device memory for configuration
	Real* dU; // to store the whole configuration (and the best copy)
	cudaMalloc( &dU, arraySize*sizeof(Real) );

	Real* dUup;
	cudaMalloc( &dUup, timesliceArraySize*sizeof(Real) ); // working copy of current timeslice t
	Real* dUdw;
	cudaMalloc( &dUdw, timesliceArraySize*sizeof(Real) ); // working copy of timeslice t-1

	// host memory for the neighbour table
	lat_index_t* nnt = (lat_index_t*)malloc( s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof(lat_index_t) );

	// device memory for the timeslice neighbour table
	lat_index_t *dNnt;
	cudaMalloc( &dNnt, s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof( lat_index_t ) );


	double* dPlaquette;
	cudaMalloc( &dPlaquette, s.getLatticeSize()*sizeof(double) );



	// initialise the timeslice neighbour table
	initNeighbourTable( nnt );
	// copy neighbour table to device
	cudaMemcpy( dNnt, nnt, s.getLatticeSizeTimeslice()*(2*(Ndim-1))*sizeof( lat_index_t ), cudaMemcpyHostToDevice );

	int threadsPerBlock = 32*8; // 32 sites are updated within a block (8 threads are needed per site)
	int numBlocks = s.getLatticeSizeTimeslice()/2/32; // // half of the lattice sites (a parity) are updated in a kernel call

	allTimer.start();

	cudaFuncSetCacheConfig( orStep, cudaFuncCachePreferL1 );
	
	// instantiate GaugeFixingStats object
//	lat_coord_t *devicePointerToSize;
//	cudaError_t error = cudaGetSymbolAddress( (void**)&devicePointerToSize, "dSize" );

	GaugeFixingStats<Ndim-1,Nc,COULOMB,AVERAGE> gaugeStats( dUup, &HOST_CONSTANTS::SIZE[1] );

//	cout << "get symbol adress error:" << cudaGetErrorString( error ) << endl;

	double totalKernelTime = 0;

	long totalStepNumber = 0;

	for( int i = fileStartnumber; i < fileStartnumber+nconf; i++ )
	{

		stringstream filename(stringstream::out);
		filename << fileBasename << setw( fileNumberformat ) << setfill( '0' ) << i << fileEnding;
//		filename << "/home/vogt/configs/STUDIENARBEIT/N32/config_n32t32beta570_sp" << setw( 4 ) << setfill( '0' ) << i << ".vogt";
		cout << "loading " << filename.str() << " as " << fileType << endl;

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
//		Real polBefore = calculatePolyakovLoopAverage( U );

		// copying configuration ...
		cudaMemcpy( dU, U, arraySize*sizeof(Real), cudaMemcpyHostToDevice );

		// calculate and print the gauge quality
		printf( "i:\t\tgff:\t\tdA:\n");
		gaugeStats.generateGaugeQuality();



		if( doMc )
		{
			for( int i = 0; i < 10000; i++ )
			{
				for( int t = 0; t < s.size[0]; t++ )
				{
					int tDw = (t==0)?s.size[0]-1:t-1;
					int tUp = (t==s.size[0]-1)?0:t+1;

					heatbathStep<<<s.getLatticeSizeTimeslice()/32/2,32>>>(&dU[tDw*timesliceArraySize], &dU[t*timesliceArraySize], &dU[tUp*timesliceArraySize], dNnt, 2.15, 0, 2*(i*HOST_CONSTANTS::SIZE[0]+t) );
					heatbathStep<<<s.getLatticeSizeTimeslice()/32/2,32>>>(&dU[tDw*timesliceArraySize], &dU[t*timesliceArraySize], &dU[tUp*timesliceArraySize], dNnt, 2.15, 1, 2*(i*HOST_CONSTANTS::SIZE[0]+t) );
				}
			}
		}






		Chronotimer kernelTimer;
		kernelTimer.reset();
		kernelTimer.start();
		for( int t = 0; t < s.size[0]; t++ )
		{
			cout << "t = " << t << endl;
			Real bestGff = 0;

			int tDw = (t==0)?s.size[0]-1:t-1;

			cudaMemcpy( dUup, &dU[t*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToDevice );
			cudaMemcpy( dUdw, &dU[tDw*timesliceArraySize], timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToDevice );

//			gaugeStats.setPointer( &dU[t*timesliceArraySize] );


			for( int copy = 0; copy < gaugeCopies; copy++ )
			{
				if( !noRandomTrafo )
				{
					randomTrafo<<<s.getLatticeSizeTimeslice()/32/2,32>>>(dUup, dUdw, dNnt, 0, 2*(t*gaugeCopies+copy) );
					randomTrafo<<<s.getLatticeSizeTimeslice()/32/2,32>>>(dUup, dUdw, dNnt, 1, 2*(t*gaugeCopies+copy)+1 );
				}

				// SIMULATED ANNEALING

				float temperature = saMax;
				float tempStep = (saMax-saMin)/(float)saSteps;
				for( int i = 0; i < saSteps; i++ )
				{

					saStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 0, temperature, i*2 );
					saStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 1, temperature, i*2+1 );

					for( int mic = 0; mic < microupdates; mic++ )
					{
						microStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 0 );
						microStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 1 );
//						gaugeStats.generateGaugeQuality();
//						printf( "%f\t\t%1.10f\t\t%e\n", temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

					}
	//				saStep<<<numBlocks,threadsPerBlock>>>(&dU[t*timesliceArraySize], &dU[tDw*timesliceArraySize], dNnt, 0, temperature, i*2 );
	//				saStep<<<numBlocks,threadsPerBlock>>>(&dU[t*timesliceArraySize], &dU[tDw*timesliceArraySize], dNnt, 1, temperature, i*2+1 );


					if( i % orCheckPrec == 0 )
					{

					// check the current gauge quality
						gaugeStats.generateGaugeQuality();
						printf( "%f\t\t%1.10f\t\t%e\n", temperature, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

			//			calculatePlaquette<<<s.getLatticeSize()/32,32>>>( dU, dNn, dPlaquette );
			//			printPlaquette<<<1,1>>>(dPlaquette );

						if( gaugeStats.getCurrentA() < orPrecision ) break;
					}
					temperature -= tempStep;
				}





				// OVERRELAXATION
				for( int i = 0; i < orMaxIter; i++ )
				{
	//				cout << "Up/Dw: " << &dU[t*timesliceArraySize] << "/" <<  &dU[tDw*timesliceArraySize] << endl;
	//				orStepSingleThread<<<s.getLatticeSizeTimeslice()/32/2,32>>>(&dU[t*timesliceArraySize], &dU[tDw*timesliceArraySize], dNnt, 0, orParameter );
	//				orStepSingleThread<<<s.getLatticeSizeTimeslice()/32/2,32>>>(&dU[t*timesliceArraySize], &dU[tDw*timesliceArraySize], dNnt, 1, orParameter );
					orStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 0, orParameter );
					orStep<<<numBlocks,threadsPerBlock>>>(dUup, dUdw, dNnt, 1, orParameter );

					if( i % orCheckPrec == 0 )
					{

					// check the current gauge quality
						gaugeStats.generateGaugeQuality();
						printf( "%d\t\t%1.10f\t\t%e\n", i, gaugeStats.getCurrentGff(), gaugeStats.getCurrentA() );

			//			calculatePlaquette<<<s.getLatticeSize()/32,32>>>( dU, dNn, dPlaquette );
			//			printPlaquette<<<1,1>>>(dPlaquette );

						if( gaugeStats.getCurrentA() < orPrecision ) break;
					}

					totalStepNumber++;
				}
				if( gaugeStats.getCurrentGff() > bestGff )
				{
					cout << "FOUND BETTER COPY" << endl;
					bestGff = gaugeStats.getCurrentGff();
					cudaMemcpy( &dU[t*timesliceArraySize], dUup, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToDevice );
					cudaMemcpy( &dU[tDw*timesliceArraySize], dUdw, timesliceArraySize*sizeof(Real), cudaMemcpyDeviceToDevice );
				}
				else
				{
					cout << "NO BETTER COPY" << endl;
				}
			}
//			exit(1);
			restoreSecondLine<<<s.getLatticeSizeTimeslice()/32,32>>>( &dU[t*timesliceArraySize], dNnt );
		}
		cudaThreadSynchronize();
		kernelTimer.stop();
		cout << "kernel time for config: " << kernelTimer.getTime() << " s"<< endl;
		totalKernelTime += kernelTimer.getTime();
		cudaMemcpy( U, dU, arraySize*sizeof(Real), cudaMemcpyDeviceToHost );

//		cout << "Polyakov loop: " << polBefore << " - " << calculatePolyakovLoopAverage( U ) << endl;



		stringstream outname(stringstream::out);
		outname << fileBasename << "gaugefixed_"<< setw( fileNumberformat ) << setfill( '0' ) << i << fileEnding;
		switch( fileType )
		{
		case VOGT:
			loadOk = lfVogt.save( s, outname.str(), U );
			break;
		case PLAIN:
			loadOk = lfPlain.save( s, outname.str(), U );
			break;
		case HEADERONLY:
			loadOk = lfHeaderOnly.save( s, outname.str(), U );
			break;
		default:
			cout << "Filetype not set to a known value. Exiting";
			exit(1);
		}
		if( loadOk ) cout << "Wrote config to file " << outname.str() << endl;
	}

	allTimer.stop();
	cout << "total time: " << allTimer.getTime() << " s" << endl;
	cout << "total kernel time: " << totalKernelTime << " s" << endl;

//	cout << (double)((long)2253*(long)s.getLatticeSize()*(long)totalStepNumber)/totalKernelTime/1.0e9 << " GFlops at "
//				<< (double)((long)192*(long)s.getLatticeSize()*(long)(totalStepNumber)*(long)sizeof(Real))/totalKernelTime/1.0e9 << "GB/s memory throughput." << endl;



}
