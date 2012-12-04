/**
 * GaugeFixingStats.hxx
 *
 *  Created on: June 13, 2012
 *      Author: schroeck
 *
 * TODO parallel reduction:
 * 	- not well adapted to smaller lattices (landau <= 24^4, coulomb <= 32^4)
 * 	- uses global memory: sould use shared memory within the block and only write the result of the local sum to global memory
 * - choose AVERAGE/MAX StoppingCrit at compile time?
 *
 */
#ifndef GAUGEFIXINGSTATS_HXX_
#define GAUGEFIXINGSTATS_HXX_

#include "GlobalConstants.h"
#include "../lattice/datatype/datatypes.h"
#include "../lattice/cuda/cuda_host_device.h"
#include "../lattice/access_pattern/StandardPattern.hxx"
#include "../lattice/access_pattern/GpuCoulombPattern.hxx"
#include "../lattice/access_pattern/GpuLandauPattern.hxx"
#include "../lattice/SiteCoord.hxx"
#include "../lattice/SiteIndex.hxx"
#include "../lattice/Link.hxx"

#include <iostream>


/**
 * TODO place somewhere else
 */
CUDA_HOST_DEVICE bool isPowerOfTwo( unsigned long x )
{
    return (x & (x - 1)) == 0;
}

__device__ inline double cuFabs( double a )
{
	return (a>0)?(a):(-a);
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> class GaugeFixingStats
{
public:
	GaugeFixingStats( Real *U, const lat_coord_t *size  );
	GaugeFixingStats( Real *U, const lat_coord_t *size, double prec  );
	~GaugeFixingStats();
	void setGaugePrecision( double  prec );
	double getGaugePrecision();
	double getCurrentGff();
	double getCurrentA();
 	void generateGaugeQuality();
 	void setPointer( Real*U );
private:
	Real *U;
	// device memory for collecting the parts of the gauge fixing functional and divA
	double *dGff;
	double *dA;
	double currentGff;
	double currentA;
	double reqPrec;
	int redBlockSize;
	lat_coord_t *dSize; // pointer to lattice size array on device constant memory
	const lat_coord_t *size; // pointer to lattice size on host memory
	SiteCoord<Ndim,FULL_SPLIT> site;

	void initReductionBlockSize();
};

template<int Ndim, int Nc, class GType, StoppingCrit ma> GaugeFixingStats<Ndim, Nc, GType, ma>::GaugeFixingStats( Real *U, const lat_coord_t *size) : site(size)
{
	this->U = U;
	this->size = size;
	this->dSize = DEVICE_CONSTANTS::SIZE;

	cudaMalloc( &dGff, site.getLatticeSize()*sizeof(double) );
	cudaMalloc( &dA,   site.getLatticeSize()*sizeof(double) );

//	std::cout << "size[0]:" << size[0] << std::endl;
//	std::cout << "size[1]:" << size[1] << std::endl;
//	std::cout << "size[2]:" << size[2] << std::endl;
//	std::cout << "LAT:" << site.getLatticeSize() << std::endl;

	redBlockSize=512;
	initReductionBlockSize();
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> GaugeFixingStats<Ndim, Nc, GType, ma>::GaugeFixingStats( Real *U, const lat_coord_t *size, double prec ) : site(size)
{
	this->U = U;
	this->size = size;
	this->dSize = DEVICE_CONSTANTS::SIZE;
	this->reqPrec = prec;
//	site(size);

	cudaMalloc( &dGff, site.getLatticeSize()*sizeof(double) );
	cudaMalloc( &dA,   site.getLatticeSize()*sizeof(double) );

	redBlockSize=512;
	initReductionBlockSize();
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> GaugeFixingStats<Ndim, Nc, GType, ma>::~GaugeFixingStats()
{
	cudaFree( &dGff );
	cudaFree( &dA );
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> void GaugeFixingStats<Ndim, Nc, GType, ma>::setPointer( Real* U )
{
	this->U = U;
	std::cout << "pointer set to " << U << std::endl;

}

/**
 * Choose reduction block size such that latticesize/redBlockSize/redBlockSize is an integer and redBlockSize is a power of 2
 *
 * For good performance we want the block size to be at least 32 (one warp).
 * TODO 24^3*T coulomb has redBlockSize = 16 => make reduction more flexible (skip reduction stages).
 */
template<int Ndim, int Nc, class GType, StoppingCrit ma> void GaugeFixingStats<Ndim, Nc, GType, ma>::initReductionBlockSize()
{
	int pot = 1;
	for( int i = 1; i < 10; i++ )
	{
		pot *= 2;
		if( site.getLatticeSize() % (pot*pot) == 0 )
		{
			redBlockSize = pot;
		}
		else
		{
			break;
		}
	}
	printf( "We chose reduction block size = %d\n", redBlockSize );
	if( isPowerOfTwo( site.getLatticeSize()/redBlockSize/redBlockSize) ) printf( "We use parallel reduction in last reduction step\n" );
	else printf( "We can't use parallel reduction in last reduction step\n" );
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> void GaugeFixingStats<Ndim, Nc, GType, ma>::setGaugePrecision( double  prec )
{
	reqPrec = prec;
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> double GaugeFixingStats<Ndim, Nc, GType, ma>::getGaugePrecision()
{
	return reqPrec;
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> double GaugeFixingStats<Ndim, Nc, GType, ma>::getCurrentGff()
{
	return currentGff;
}

template<int Ndim, int Nc, class GType, StoppingCrit ma> double GaugeFixingStats<Ndim, Nc, GType, ma>::getCurrentA()
{
	return currentA;
}



//template<int Ndim, int Nc, GaugeType lc, StoppingCrit ma>  __global__ void generateGaugeQualityPerSite( Real *U, double *dGff, double *dA, lat_coord_t *size );
__global__ void reduce1GaugeQuality( double *dGff, double *dA, StoppingCrit ma );
__global__ void reduce2GaugeQuality( double *dGff, double *dA, StoppingCrit ma );
__global__ void reduce3GaugeQuality( double *dGff, double *dA, int redBlockSize, StoppingCrit ma );
__global__ void reduceGaugeQuality( double *dGff, double *dA, int latticeSize, StoppingCrit ma );
__device__ double average_or_max( double a, double b, StoppingCrit ma);



template<int Ndim, int Nc, class GType, StoppingCrit ma> void GaugeFixingStats<Ndim,Nc,GType,ma>::generateGaugeQuality()
{

	GType::generateGaugeQualityPerSite(site.getLatticeSize()/32, 32, U, dGff, dA );
//	generateGaugeQualityPerSite<Ndim,Nc,lc,ma><<<site.getLatticeSize()/32,32>>>(U, dGff, dA, dSize);

	reduce1GaugeQuality<<<site.getLatticeSize()/redBlockSize,redBlockSize>>>(dGff, dA, ma);
	reduce2GaugeQuality<<<site.getLatticeSize()/redBlockSize/redBlockSize,redBlockSize>>>(dGff, dA, ma);
	reduce3GaugeQuality<<<1,site.getLatticeSize()/redBlockSize/redBlockSize>>>(dGff, dA, redBlockSize, ma);

//		reduceGaugeQuality<<<1,1>>>(dGff,dA,site.getLatticeSize(), ma);

	cudaMemcpy( &currentGff, dGff, sizeof(double), cudaMemcpyDeviceToHost );
	cudaMemcpy( &currentA,   dA,   sizeof(double), cudaMemcpyDeviceToHost );
	
	currentGff = currentGff*GType::getGaugeQualityPrefactorGff()/double(site.getLatticeSize());
	currentA   = currentA*GType::getGaugeQualityPrefactorA();
	if( ma == AVERAGE )
		currentA   = currentA/double(site.getLatticeSize());


	
//	if( lc == LANDAU || lc == U1xU1 )
//	{
//		currentGff = currentGff/double(site.getLatticeSize())/(double)Ndim/(double)Nc;
//		currentA   = currentA/(double)Nc;
//		if( ma == AVERAGE )
//			currentA   = currentA/double(site.getLatticeSize());
//	}
//	else if( lc == COULOMB )
//	{
//		currentGff = currentGff/double(site.getLatticeSize())/(double)(Ndim)/(double)Nc;
//		currentA   = currentA/(double)Nc;
//		if( ma == AVERAGE )
//			currentA   = currentA/double(site.getLatticeSize());
//	}
//	else if( lc == MAG )
//	{
//		currentGff = currentGff*2.0/double(site.getLatticeSize())/(double)Nc;
//		currentA   = currentA/(double)Nc/(double)Ndim/2.0;
//		if( ma == AVERAGE )
//			currentA   = currentA/double(site.getLatticeSize());
//	}
	
}


//template<int Ndim, int Nc, GaugeType lc, StoppingCrit ma> __global__ void  generateGaugeQualityPerSite( Real *U, double *dGff, double *dA, lat_coord_t *size )
//{
//	// TODO we have to copy the full code for LANDAU and COULOMB, because of different patterns. Thats why I choose to have different applications for Landau and Coulomb gauge fixing...
//
////	if(blockIdx.x * blockDim.x + threadIdx.x == 0 ) printf( "site %f\n", U[0] );
//
//	if( Nc == 2 )
//	{
//		if( lc == LANDAU )
//		{
//			typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu; // TODO check parity
//			typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
//
//			SiteCoord<Ndim,true> s(size);
//			int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//			double result = 0;
//
//			double prec1 = 0;
//			double prec2 = 0;
//			double prec3 = 0;
//
////			SU2<Quaternion<Real> > temp(locTemp);
//			for( int mu = 0; mu < 4; mu++ )
//			{
//				s.setLatticeIndex( site );
//				TLink linkUp( U, s, mu );
//				SU2<TLink> globUp( linkUp );
////				temp.assignQuaternion( globUp );
////				result += locTemp[0];
////				result += temp.get(0,0).x;
//				complex temp = globUp.get(0,0);
//				result += temp.x*2.;
//
//				prec1 += temp.y;
//				temp = globUp.get(0,1);
//				prec2 += temp.x;
//				prec3 += temp.y;
//
//				s.setNeighbour(mu,-1);
//				TLink linkDw( U, s, mu );
//				SU2<TLink> globDw( linkDw );
//
//				prec1 -= globDw.get(0,0).y;
//				temp = globDw.get(0,1);
//				prec2 -= temp.x;
//				prec3 -= temp.y;
//			}
//
//			dA[site] = cuFabs( prec1 ) + cuFabs( prec2 ) + cuFabs( prec3 );
//			dGff[site] = result;
//		}
//
//		if( lc == COULOMB )
//		{
//			typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
//			typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
//			int site = blockIdx.x * blockDim.x + threadIdx.x;
//
////			if( site == 0 ) printf( "SIZE[0]:%d\n", DEVICE_CONSTANTS::SIZE[0] );
////			if( site == 0 ) printf( "SIZE[1]:%d\n", DEVICE_CONSTANTS::SIZE[1] );
////			if( site == 0 ) printf( "SIZE[2]:%d\n", DEVICE_CONSTANTS::SIZE[2] );
////			if( site == 0 ) printf( "SIZE[3]:%d\n", DEVICE_CONSTANTS::SIZE[3] );
//
//
////			lat_coord_t tempSize[Ndim];
////			for( int i = 0; i < Ndim; i++ )
////			{
////				tempSize[i] = DEVICE_CONSTANTS::SIZE[i+1];
////			}
////
////			SiteCoord<Ndim,true> s(tempSize);
//
//
//			SiteCoord<Ndim,true> s(&DEVICE_CONSTANTS::SIZE[1]);
//
//
////			if( site == 0 ) printf( "LATSIZE:%d\n", s.getLatticeSize() );
//
//			double result = 0;
//
//			double prec1 = 0;
//			double prec2 = 0;
//			double prec3 = 0;
//
////			SU2<Quaternion<Real> > temp(locTemp);
//			for( int mu = 1; mu < 4; mu++ )
//			{
//				s.setLatticeIndex( site );
//				TLink linkUp( U, s, mu );
//				SU2<TLink> globUp( linkUp );
////				temp.assignQuaternion( globUp );
////				result += locTemp[0];
////				result += temp.get(0,0).x;
//				complex temp = globUp.get(0,0);
//				result += temp.x*2.;
//
//				prec1 += temp.y;
//				temp = globUp.get(0,1);
//				prec2 += temp.x;
//				prec3 += temp.y;
//
//				s.setNeighbour(mu-1,-1);
//				TLink linkDw( U, s, mu );
//				SU2<TLink> globDw( linkDw );
//
//				prec1 -= globDw.get(0,0).y;
//				temp = globDw.get(0,1);
//				prec2 -= temp.x;
//				prec3 -= temp.y;
//			}
//
//			dA[site] = cuFabs( prec1 ) + cuFabs( prec2 ) + cuFabs( prec3 );
//			dGff[site] = result;
//
////			dA[site] = 1;
////			dGff[site] = 1;
//
//
//		}
//	}
//
//	else if( Nc == 3 )
//	{
//		if( lc == LANDAU || lc == U1xU1 )
//		{
//			typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
//			typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
//
//			SiteCoord<Ndim,true> s(DEVICE_CONSTANTS::SIZE);
//			int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//			Matrix<complex,Nc> locMatSum;
//			SU3<Matrix<complex,Nc> > Sum(locMatSum);
//
//			Sum.zero();
//
//			// TODO check if there is a faster way to compute DELTA
//			for( int mu = 0; mu < 4; mu++ )
//			{
//				s.setLatticeIndex( site );
//
//				Matrix<complex,Nc> locMat;
//				SU3<Matrix<complex,Nc> > temp(locMat);
//
//				TLink linkUp( U, s, mu );
//				SU3<TLink> globUp( linkUp );
//
//				temp.assignWithoutThirdLine( globUp );
////				temp.projectSU3withoutThirdRow();// TODO project here?
//				globUp.assignWithoutThirdLine( temp ); // TODO
//				temp.reconstructThirdLine();
//				Sum += temp;
//
//				s.setNeighbour(mu,-1);
//				TLink linkDw( U, s, mu );
//				SU3<TLink> globDw( linkDw );
//				temp.assignWithoutThirdLine( globDw );
////				temp.projectSU3withoutThirdRow(); // TODO project here?
//				globDw.assignWithoutThirdLine( temp ); // TODO
//				temp.reconstructThirdLine();
//				Sum -= temp;
//			}
//
//			Sum -= Sum.trace()/Real(3.);
//
//			Matrix<complex,Nc> locMatSumHerm;
//			SU3<Matrix<complex,Nc> > SumHerm(locMatSumHerm);
//			SumHerm = Sum;
//			SumHerm.hermitian();
//
//			Sum -= SumHerm;
//
//			double prec = 0;
//
//			if( lc == LANDAU )
//			{
//				for( int i = 0; i < 3; i++ )
//				{
//					for( int j = 0; j < 3; j++ )
//					{
//						prec += Sum.get(i,j).abs_squared();
//					}
//				}
//			}
//			else if( lc == U1xU1 ) // only diagonal part
//			{
//				for( int j = 0; j < 3; j++ )
//				{
//					prec += Sum.get(j,j).abs_squared();
//				}
//			}
//
//			dA[site] = prec;
//
//
//			s.setLatticeIndex( site );
//			double result = 0;
//
//			Matrix<complex,Nc> locTemp;
//			SU3<Matrix<complex,Nc> > temp(locTemp);
//			for( int mu = 0; mu < 4; mu++ )
//			{
//				TLink linkUp( U, s, mu );
//				SU3<TLink> globUp( linkUp );
//				temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
//				temp.reconstructThirdLine();
//				result += temp.trace().x;
//			}
//
//			dGff[site] = result;
//		}
//		else if( lc == COULOMB )
//		{
//			// This uses a Coulomb Memory pattern (see redefinition of Size to use GpuLandauPattern in Ndim-1 dimensions
//			typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
//			typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
//
////			lat_coord_t tempSize[Ndim];
////			for( int i = 0; i < Ndim; i++ )
////			{
////				tempSize[i] = size[i+1];
////			}
//
////			SiteCoord<Ndim,true> s(tempSize);
//
//			SiteCoord<Ndim,true> s(&DEVICE_CONSTANTS::SIZE[1]);
//
//			int site = blockIdx.x * blockDim.x + threadIdx.x;
//
//			Matrix<complex,Nc> locMatSum;
//			SU3<Matrix<complex,Nc> > Sum(locMatSum);
//
//			Sum.zero();
//
//			// TODO check if there is a faster way to compute DELTA
//			for( int mu = 1; mu < 4; mu++ )
//			{
//				s.setLatticeIndex( site );
//
//				Matrix<complex,Nc> locMat;
//				SU3<Matrix<complex,Nc> > temp(locMat);
//
//				TLink linkUp( U, s, mu );
//				SU3<TLink> globUp( linkUp );
//
//				temp.assignWithoutThirdLine( globUp );
//				temp.projectSU3withoutThirdRow();// TODO project here?
//				globUp.assignWithoutThirdLine( temp ); // TODO
//				temp.reconstructThirdLine();
//				Sum += temp;
//
//				s.setNeighbour(mu-1,-1);
//				TLink linkDw( U, s, mu );
//				SU3<TLink> globDw( linkDw );
//				temp.assignWithoutThirdLine( globDw );
//				temp.projectSU3withoutThirdRow(); // TODO project here?
//				globDw.assignWithoutThirdLine( temp ); // TODO
//				temp.reconstructThirdLine();
//				Sum -= temp;
//			}
//
//			Sum -= Sum.trace()/Real(3.);
//
//			Matrix<complex,Nc> locMatSumHerm;
//			SU3<Matrix<complex,Nc> > SumHerm(locMatSumHerm);
//			SumHerm = Sum;
//			SumHerm.hermitian();
//
//			Sum -= SumHerm;
//
//			double prec = 0;
//			for( int i = 0; i < 3; i++ )
//			{
//				for( int j = 0; j < 3; j++ )
//				{
//					prec += Sum.get(i,j).abs_squared();
//				}
//			}
//
//			dA[site] = prec;
//
//
//			s.setLatticeIndex( site );
//			double result = 0;
//
//			Matrix<complex,Nc> locTemp;
//			SU3<Matrix<complex,Nc> > temp(locTemp);
//			for( int mu = 1; mu < 4; mu++ )
//			{
//				TLink linkUp( U, s, mu );
//				SU3<TLink> globUp( linkUp );
//				temp.assignWithoutThirdLine( globUp );  // TODO put this in the loop of dA to reuse globUp
//				temp.reconstructThirdLine();
//				result += temp.trace().x;
//			}
//
//			dGff[site] = result;
//		}
//	}
//	else if( lc == MAG )
//	{
//		typedef GpuLandauPattern< SiteCoord<Ndim,true>,Ndim,Nc> Gpu;
//		typedef Link<Gpu,SiteCoord<Ndim,true>,Ndim,Nc> TLink;
//
//		SiteCoord<Ndim,true> s(DEVICE_CONSTANTS::SIZE);
//		int site = blockIdx.x * blockDim.x + threadIdx.x;
//		Quaternion<Real> u,v;
//		Complex<double> X[3];
//
//		dGff[site]=0.0;
//
//		for( int mu = 0; mu < 4; mu++ )
//		{
//			s.setLatticeIndex( site );
//			TLink linkUp( U, s, mu );
//			SU3<TLink> globUp( linkUp );
////			globUp.projectSU3withoutThirdRow();// TODO project here?
//
//			s.setNeighbour(mu,-1);
//			TLink linkDw( U, s, mu );
//			SU3<TLink> globDw( linkDw );
////			globDw.projectSU3withoutThirdRow();// TODO project here?
//
//			// action:
//			for( int j=0; j<3; j++ )
//				dGff[site] += linkUp.get(j,j).abs_squared();
//
//			// precision:
//			for( int i=0; i<2; i++ )
//				for( int j=i+1; j<3; j++ )
//				{
//					u = globUp.getSubgroupQuaternion( i, j );
//					u.projectSU2();
//
//					v = globDw.getSubgroupQuaternion( i, j );
//					v.projectSU2();
//
//					X[i+j-1].x += -u[0]*u[2] + u[1]*u[3] + v[0]*v[2] + v[1]*v[3];
//					X[i+j-1].y += -u[0]*u[1] - u[2]*u[3] + v[0]*v[1] - v[2]*v[3];
//				}
//		}
//		dA[site] = X[0].abs_squared()+X[1].abs_squared()+X[2].abs_squared();
//	}
//}

__global__ void reduce1GaugeQuality( double *dGff, double *dA, StoppingCrit ma )
{
	int site = blockIdx.x*blockDim.x + threadIdx.x;
	
// 	for( int i = 0; i < latticeSize; i++ )
// 	{
// 		gff+= dGff[i];
// 		if( cuFabs(dA[i]) > temp ) temp = cuFabs(dA[i]);
// 	}
	// Perform tree-like reduction of accumulators' results.
	// blockDim.x has to be power of two at this stage 
	
	// reduce within each thread block
	for( int stride=blockDim.x/2; stride>0; stride>>=1 )
	{
		__syncthreads();
		if( threadIdx.x<stride ) 
		{
			dGff[site] += dGff[site+stride];
			  dA[site] = average_or_max( dA[site], dA[site+stride], ma );
		}
	}
}

__global__ void reduce2GaugeQuality( double *dGff, double *dA, StoppingCrit ma )
{
	int i = blockIdx.x*blockDim.x*blockDim.x + threadIdx.x*blockDim.x;
	
	for( int stride=blockDim.x/2; stride>0; stride>>=1 )
	{
		__syncthreads();
		if( threadIdx.x<stride ) 
		{
			dGff[i] += dGff[i+stride*blockDim.x];
			  dA[i] = average_or_max( dA[i], dA[i+stride*blockDim.x], ma );
		}
	}
// 	for( int i=bs; i<latticeSize; i+=bs )
// 	{
// 		dGff[0] += dGff[i];
// 		dA[0]   += dA[i];
// 	}
}


__global__ void reduce3GaugeQuality( double *dGff, double *dA, int redBlockSize, StoppingCrit ma )
{
	int bs = redBlockSize*redBlockSize;
	int i = bs*threadIdx.x;

	if( isPowerOfTwo( blockDim.x ) ) // we can use parallel reduction
	{
		for( int stride=blockDim.x/2; stride>0; stride>>=1 )
		{
			__syncthreads();
			if( threadIdx.x<stride )
			{
				dGff[i] += dGff[i+stride*bs];
				  dA[i] = average_or_max( dA[i], dA[i+stride*bs], ma );
			}
		}
	}
	else // serialise the last iteration (only one thread is active)
	{
		if (i == 0 )
			for( int j = 1; j < blockDim.x; j ++ )
			{
				dGff[0] += dGff[0+bs*j];
				  dA[0] = average_or_max( dA[0], dA[0+bs*j], ma );
			}
	}
}


__global__ void reduceGaugeQuality( double *dGff, double *dA, int latticeSize, StoppingCrit ma )
{
	double gff = 0;
	double A = 0;
	for( int i = 0; i < latticeSize; i++ )
	{
		gff+= dGff[i];
		A = average_or_max( A, dA[i], ma );
//		if( cuFabs(dA[i]) > temp ) temp = cuFabs(dA[i]);
	}

	dGff[0] = gff;
	dA[0] = A;
}

__device__ double average_or_max( double a, double b, StoppingCrit ma)
{
	if( ma == AVERAGE ) return a+b;
	else if ( ma == MAX ) return (a>b?a:b);
	else return 0;
}



#endif /* GAUGEFIXINGSTATS_HXX_ */
