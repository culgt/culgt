/**
 * rootm.h
 *
 *  Created on: Jun 30, 2014
 *      Author: vogt
 */

#ifndef ROOTM_H_
#define ROOTM_H_

#include "itpp/itbase.h"
#include <complex>
#include "../lattice/parameterization_types/SU2Vector4.h"
#include "../lattice/LocalLink.h"
#include "../lattice/parameterization_types/SUNRealFull.h"

using culgt::LocalLink;
using culgt::SUNRealFull;

typedef std::complex<double> complex_double;



class Rootm
{
public:
	static itpp::cmat root( itpp::cmat A, int p );
	template<typename LinkType> static LinkType root( LinkType Q, int p )
	{
		LocalLink<SUNRealFull<LinkType::PARAMTYPE::NC,typename LinkType::PARAMTYPE::REALTYPE> > link;
		link = Q;

		itpp::cmat A(LinkType::PARAMTYPE::NC,LinkType::PARAMTYPE::NC);

		for( int i = 0; i < LinkType::PARAMTYPE::NC; i++ )
		{
			for( int j = 0; j < LinkType::PARAMTYPE::NC; j++ )
			{
				A(i,j) = complex_double( link[2*(i*LinkType::PARAMTYPE::NC+j)], link[2*(i*LinkType::PARAMTYPE::NC+j)+1] );
			}
		}

		A = root( A, p );


		for( int i = 0; i < LinkType::PARAMTYPE::NC; i++ )
		{
			for( int j = 0; j < LinkType::PARAMTYPE::NC; j++ )
			{
				link[2*(i*LinkType::PARAMTYPE::NC+j)] = std::real( A(i,j) );
				link[2*(i*LinkType::PARAMTYPE::NC+j)+1] = std::imag( A(i,j) );
			}
		}

		Q = link;
		return Q;
	}
};

itpp::cmat Rootm::root( itpp::cmat A, int p )
{
	// assume A is square matrix
	int n = A.cols();

	itpp::cmat Q(n, n);
	itpp::cmat T(n, n);
	itpp::cmat U(n, n);
	itpp::cmat R(n,(p-2)*n);


	itpp::schur( A, Q, T );

	U.zeros();
	R.zeros();

	for( int i = 1; i <= n; i++ )
	{
		U(i-1,i-1) = pow( T(i-1,i-1), 1/double(p) );
		for( int a = 1; a <= p-2; a++ )
		{
			R(i-1,(a-1)*n+i-1) = pow(U(i-1,i-1),(a+2));
		}
	}


	complex_double sum1;
	complex_double sum2;
	complex_double sum3;
	complex_double sum4;
	complex_double sum5;
	complex_double sum6;
	complex_double sum7;

	complex_double sum;
	for( int c = 1; c <= n-1; c++ ) // using matlab notation and subtract 1 when accesing matrix elements
	{
		for( int i = 1; i <= n-c; i++ )
		{
			sum1 = 0;
			for( int d = 1; d <= p-2; d++ )
			{
				sum2 = 0;
				for( int k = i+1; k <= i+c-1; k++ )
				{
					sum2 += U(i-1,k-1)*R(k-1,(d-1)*n+i+c-1);
				}
				sum1 += pow( U(i-1,i-1), (p-2-d))*sum2;
			}
			sum3 = 0;
			for( int k = i+1; k <= i+c-1; k++ )
			{
				sum3 += U(i-1,k-1)*U(k-1,i+c-1);
			}
			sum1 += pow(U(i-1,i-1), (p-2))*sum3;

			sum4 = 0;
			for( int j = 0; j <= p-1; j++ )
			{
				sum4 += pow( U(i-1,i-1), j )*pow(U(i+c-1,i+c-1), p-1-j);
			}

			U( i-1, i+c-1 ) = ( T(i-1,i+c-1) - sum1)/sum4;

			for( int q = 1; q <= p-2; q++ )
			{
				sum5 = 0;
				for( int g = 0; g <= q; g++ )
				{
					sum5+= pow(U(i-1,i-1), g) * pow( U(i+c-1,i+c-1), q-g);
				}
				sum6 = 0;
				for( int h = 1; h <= q-1; h++ )
				{
					sum7 = 0;
					for( int w = i+1; w <= i+c-1; w++ )
					{
						sum7 += U(i-1,w-1)*R(w-1,(h-1)*n+i+c-1);
					}
					sum6 += pow(U(i-1,i-1), q-1-h)*sum7;
				}
				sum = sum6 + pow(U(i-1,i-1),q-1)*sum3;

				R(i-1,(q-1)*n+i+c-1 ) = U(i-1, i+c-1)*sum5 + sum;
			}
		}
	}

	return Q*U*hermitian_transpose(Q);
}





#endif /* ROOTM_H_ */
