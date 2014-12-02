/**
 * Cardano3x3.h
 *
 * T can be float, double, Complex<float>, Complex<double>
 *
 *  Created on: Jul 4, 2014
 *      Author: vogt
 */

#ifndef CARDANO3X3_H_
#define CARDANO3X3_H_
#include "cudacommon/cuda_host_device.h"
#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include <iostream>
#include <assert.h>
#include <stdio.h>

namespace culgt
{

template<typename T> class Cardano3x3
{
public:
	CUDA_HOST_DEVICE inline Cardano3x3( Matrix<T,3>& mat ) : mat(mat)
	{
		compute();
	}

	CUDA_HOST_DEVICE inline T getEigenvalue( int i )
	{
		return p*beta[i]+q;
	}

	CUDA_HOST_DEVICE inline Vector<T,3> getEigenvector( int i )
	{
		T b[2];
		int counter = 0;
		for( int j = 0; j < 3; j++ )
		{
			if( i != j )
			{
				b[counter] = beta[j];
				counter++;
			}
		}

		Matrix<T,3> tmp = (B-b[0]);
		tmp *= (B-b[1]);
		Vector<T,3> eigvec = tmp.col(0);

//		if( isnan(  (eigvec / eigvec.norm())(0).x ) || isnan(  (eigvec / eigvec.norm())(0).y ) )
//		{
//			printf( "eigvec %f+i*%f \t%f+i*%f \t%f+i*%f \t \n", eigvec(0).x,eigvec(0).y,eigvec(1).x,eigvec(1).y,eigvec(2).x,eigvec(2).y);
//			printf( "b[0] %f+i*%f\n", b[0].x, b[0].y);
//			printf( "b[1] %f+i*%f\n", b[1].x, b[1].y);
//		}

		return eigvec / eigvec.norm();
	}

	CUDA_HOST_DEVICE inline Matrix<T,3> getEigenvectors()
	{
		Matrix<T,3> result;
		for( int i = 0; i < 3; i++ )
		{
			result.setCol( i, getEigenvector( i ) );
		}
		return result;
	}

	CUDA_HOST_DEVICE inline void logAssumeUnitary()
	{
		Matrix<T,3> S;
		S = getEigenvectors();
		Matrix<T,3> Sinv = S;
		Sinv.hermitian();

		Matrix<T,3> temp;
		temp = Sinv*mat;
		temp*=S;

		for( int i = 0; i < 3; i++ )
		{
//			if( isnan(  log(temp(i,i)).x ) || isnan(  log(temp(i,i)).y ) )
//			{
//				printf("%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \n", Sinv(0,0).x, Sinv(0,0).y, Sinv(0,1).x, Sinv(0,1).y, Sinv(0,2).x, Sinv(0,2).y, Sinv(1,0).x, Sinv(1,0).y, Sinv(1,1).x, Sinv(1,1).y, Sinv(1,2).x, Sinv(1,2).y, Sinv(2,0).x, Sinv(2,0).y, Sinv(2,1).x, Sinv(2,1).y, Sinv(2,2).x, Sinv(2,2).y );
//				printf("%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \n", S(0,0).x, S(0,0).y, S(0,1).x, S(0,1).y, S(0,2).x, S(0,2).y, S(1,0).x, S(1,0).y, S(1,1).x, S(1,1).y, S(1,2).x, S(1,2).y, S(2,0).x, S(2,0).y, S(2,1).x, S(2,1).y, S(2,2).x, S(2,2).y );
//				printf( "logDiag.x is %f\n", log(temp(i,i)).x);
//				printf( "logDiag.y is %f\n", log(temp(i,i)).y);
//				printf( "diag.x is %f\n", temp(i,i).x);
//				printf( "diag.y is %f\n", temp(i,i).y);
//
//			}
			temp(i,i) = log( temp(i,i) );
		}

		mat = S*temp;
		mat *= Sinv;
	}

private:
	Matrix<T,3>& mat;
	Matrix<T,3> B;
	T beta[3]; // eigenvalues
	T p,q;

	/**
	 * computes the eigenvalues
	 * @return
	 */
	CUDA_HOST_DEVICE inline void compute()
	{
		T off = 0.;
		for( int i = 0; i < 3; i++ )
			for( int j = 0; j < 3; j++ )
			{
				if( i != j ) off += mat(i,j);
			}

		// we are already diagonal -> save eigenvalues
		if( fabs(off) == 0 )
		{
			p = 1.;
			q = 0.;
			for( int i = 0; i < 3; i++ )
			{
				beta[i] = mat(i,i);
			}
		}
		else
		{
			q = mat.trace()/3.;
			Matrix<T,3> tmp;
			tmp = mat-q;
			tmp *= tmp;
			p = sqrt( tmp.trace()/6. );

//			if( isnan(  p.x ) ) printf( "p is nan\n");
//			if( p.x == 0 )
//			{
//				printf("%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \n", mat(0,0).x, mat(0,0).y, mat(0,1).x, mat(0,1).y, mat(0,2).x, mat(0,2).y, mat(1,0).x, mat(1,0).y, mat(1,1).x, mat(1,1).y, mat(1,2).x, mat(1,2).y, mat(2,0).x, mat(2,0).y, mat(2,1).x, mat(2,1).y, mat(2,2).x, mat(2,2).y );
//
////				printf( "p is 0\n");
//
//				// TODO how to handle p = 0???
//			}
//			else
//			{
				B = mat-q;
				B /= p;
//			}


//#ifndef __CUDA_ARCH__
//			if( fabs( B.det() ) > 2. )
//			{
//				std::cout << "algorithm not implemented for det > 2" << std::endl;
//				assert(false);
//			}
//#endif
			T r = B.det()/2.;

//			if( isnan(  fabs( B.det() ) ) ) printf( "det is nan\n");

			T acosr = acos(r);
//			if( isnan(  fabs( acosr ) ) ) printf( "acosr is nan\n");

			beta[0] = 2.*cos( 1./3.*acosr );
			beta[1] = 2.*cos( 1./3.*acosr + 2.*::acos(-1.)/3.);
			beta[2] = 2.*cos( 1./3.*acosr + 4.*::acos(-1.)/3. );
//			if( isnan(  fabs( beta[0] ) )  || isinf(fabs( beta[0] ) ) )
//			{
//				printf("%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \n", mat(0,0).x, mat(0,0).y, mat(0,1).x, mat(0,1).y, mat(0,2).x, mat(0,2).y, mat(1,0).x, mat(1,0).y, mat(1,1).x, mat(1,1).y, mat(1,2).x, mat(1,2).y, mat(2,0).x, mat(2,0).y, mat(2,1).x, mat(2,1).y, mat(2,2).x, mat(2,2).y );
//				printf("%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \t\n%f+i*%f \t%2.15f+i*%2.15f \t%2.15f+i*%2.15f \n", B(0,0).x, B(0,0).y, B(0,1).x, B(0,1).y, B(0,2).x, B(0,2).y, B(1,0).x, B(1,0).y, B(1,1).x, B(1,1).y, B(1,2).x, B(1,2).y, B(2,0).x, B(2,0).y, B(2,1).x, B(2,1).y, B(2,2).x, B(2,2).y );
//				printf( "r is %2.15f+i*%2.15f\n", r.x, r.y);
//				printf( "acosr is %2.15f+i*%2.15f\n", acosr.x, acosr.y);
//				T check = r + T::I()*sqrt( T(1.,0.) - r*r );
//				T check2 = sqrt( T(1.,0.) - r*r );
//				T check3 = T::I()*sqrt( T(1.,0.) - r*r );
//				printf( "r + i*sqrt( 1. - r*r ) is %2.15f+i*%2.15f\n", check.x, check.y);
//				printf( "sqrt( 1. - r*r ) is %2.15f+i*%2.15f\n", check2.x, check2.y);
//				printf( "i*sqrt( 1. - r*r ) is %2.15f+i*%2.15f\n", check3.x, check3.y);
//				check3 += r;
//				printf( "i*sqrt( 1. - r*r )+r is %2.15f+i*%2.15f\n", check3.x, check3.y);
//				printf( "r is %2.15f+i*%2.15f\n", r.x, r.y);
//				T check4 = T::I()*r + sqrt( T(1.,0.) - r*r );
//				printf( "i*r+sqrt( 1. - r*r ) is %2.15f+i*%2.15f\n", check4.x, check4.y);
//				printf( "log(i*r+sqrt( 1. - r*r )) is %2.15f+i*%2.15f\n", log(check4).x, log(check4).y);
//				check4 = log(check4)*T::I()*-1.;
//				printf( "-i*log(i*r+sqrt( 1. - r*r )) is %2.15f+i*%2.15f\n",  check4.x, check4.y);
//
//				printf( "p is %2.15f+i*%2.15f\n", p.x, p.y);
//				printf( "q is %2.15f+i*%2.15f\n", q.x, q.y);
//			}
//			if( isnan(  fabs( beta[1] ) ) ) printf( "beta1 is nan\n");
//			if( isnan(  fabs( beta[2] ) ) ) printf( "beta2 is nan\n");
		}
	}
};
}

#endif /* CARDANO3X3_H_ */
