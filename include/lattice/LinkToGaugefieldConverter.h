/**
 * LinkToGaugefieldConvert.h
 *
 *  Created on: Jun 18, 2014
 *      Author: vogt
 */

#ifndef LINKTOGAUGEFIELDCONVERTER_H_
#define LINKTOGAUGEFIELDCONVERTER_H_

#include "../cudacommon/cuda_host_device.h"
#include <assert.h>
#include "parameterization_types/SU2Vector4.h"
#include "parameterization_types/SU2Real4.h"
#include "parameterization_types/SUNComplexFull.h"
#include "LocalLink.h"
#include "../util/math/GellMannMatrices.h"
#include "eigensolver/Matrix.h"
#include "eigensolver/Cardano3x3.h"
#include <boost/algorithm/string.hpp>
#include "gaugefieldtypes.h"

using boost::algorithm::iequals;

namespace culgt
{

template<int Nc, gaugefieldtype::GaugefieldType GFType> class LinkToGaugefieldConverter
{
};

template<> class LinkToGaugefieldConverter<2, gaugefieldtype::LINEAR>
{
public:
	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( T* A, LinkType& link )
	{
		LocalLink<SU2Vector4<typename LinkType::PARAMTYPE::REALTYPE> > linkConverted;
		linkConverted = link;

		A[0] = 2.*(T)linkConverted.get(0).w;
		A[1] = 2.*(T)linkConverted.get(0).z;
		A[2] = 2.*(T)linkConverted.get(0).y;
	}

	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( LinkType& link,  T* A )
	{
		LocalLink<SU2Vector4<typename LinkType::PARAMTYPE::REALTYPE> > linkConverted = link;

		typename Real4<typename LinkType::PARAMTYPE::REALTYPE>::VECTORTYPE temp;

		temp.w = A[0]*.5;
		temp.z = A[1]*.5;
		temp.y = A[2]*.5;
		temp.x = ::sqrt( 1. - ( temp.y*temp.y + temp.z*temp.z + temp.w*temp.w ) );

		linkConverted.set( 0, temp );
		link = linkConverted;
	}
};

template<> class LinkToGaugefieldConverter<2, gaugefieldtype::LOGARITHMIC>
{
public:
	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( T* A, LinkType& link )
	{
		LocalLink<SU2Vector4<typename LinkType::PARAMTYPE::REALTYPE> > linkConverted;
		linkConverted = link;

		T rUAbs = rsqrt( linkConverted.get(0).w*linkConverted.get(0).w + linkConverted.get(0).z*linkConverted.get(0).z + linkConverted.get(0).y*linkConverted.get(0).y );
		T acosU0 = ::acos( linkConverted.get(0).x );

		A[0] = 2.*(T)linkConverted.get(0).w*rUAbs*acosU0;
		A[1] = 2.*(T)linkConverted.get(0).z*rUAbs*acosU0;
		A[2] = 2.*(T)linkConverted.get(0).y*rUAbs*acosU0;
	}

	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( LinkType& link, T* A )
	{
		T aAbs = ::sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );

		LocalLink<SU2Vector4<typename LinkType::PARAMTYPE::REALTYPE> > linkConverted;
		typename Real4<typename LinkType::PARAMTYPE::REALTYPE>::VECTORTYPE temp;

		temp.x = ::cos( .5 * aAbs );
		temp.w = ::sin( .5 * aAbs ) * A[0]/aAbs;
		temp.z = ::sin( .5 * aAbs ) * A[1]/aAbs;
		temp.y = ::sin( .5 * aAbs ) * A[2]/aAbs;

		linkConverted.set(0, temp );
		link = linkConverted;
	}
};


template<> class LinkToGaugefieldConverter<3, gaugefieldtype::LINEAR>
{
public:
	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( T* A, LinkType& link )
	{
		typedef typename LinkType::PARAMTYPE::REALTYPE REALTYPE;
		LocalLink<SUNComplexFull<3,REALTYPE> > linkConverted;
		LocalLink<SUNComplexFull<3,REALTYPE> > linkConvertedHerm;
		linkConverted = link;
		linkConvertedHerm = linkConverted;
		linkConvertedHerm.hermitian();

		linkConverted -= linkConvertedHerm;

		Complex<T> trace = linkConverted.trace();
		linkConvertedHerm.identity();
		linkConvertedHerm *= (trace/3.);

		linkConverted -= linkConvertedHerm;
		linkConverted *= 1./(2.*Complex<T>::I());

		for( int i = 0; i < 8; i++ )
		{
			linkConvertedHerm = linkConverted;
			linkConvertedHerm *= GellMannMatrices<LocalLink<SUNComplexFull<3,REALTYPE> > >::get(i+1);
			A[i] = linkConvertedHerm.reTrace();
		}
	}
};

template<> class LinkToGaugefieldConverter<3, gaugefieldtype::LOGARITHMIC>
{
public:
	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( T* A, LinkType& link )
	{
		typedef typename LinkType::PARAMTYPE::REALTYPE REALTYPE;
		LocalLink<SUNComplexFull<3,REALTYPE> > linkConverted;
		LocalLink<SUNComplexFull<3,REALTYPE> > temp;
		linkConverted = link;

		Matrix<Complex<T>,3>* mat = reinterpret_cast<Matrix<Complex<T>,3>*>( &linkConverted );

		Cardano3x3<Complex<T> > eigensolver( *mat );
		eigensolver.logAssumeUnitary();

		*mat*=Complex<T>(0,-1);

		for( int i = 0; i < 8; i++ )
		{
			temp = linkConverted;
			temp *= GellMannMatrices<LocalLink<SUNComplexFull<3,REALTYPE> > >::get(i+1);
			A[i] = temp.reTrace();
		}
	}

//	template<typename LinkType, typename T> static CUDA_HOST_DEVICE void convert( LinkType& link,  T* A )
//	{
//		LocalLink<SU2Vector4<typename LinkType::PARAMTYPE::REALTYPE> > linkConverted = link;
//
//		typename Real4<typename LinkType::PARAMTYPE::REALTYPE>::VECTORTYPE temp;
//
//		temp.w = A[0]*.5;
//		temp.z = A[1]*.5;
//		temp.y = A[2]*.5;
//		temp.x = sqrt( 1. - ( temp.y*temp.y + temp.z*temp.z + temp.w*temp.w ) );
//
//		linkConverted.set( 0, temp );
//		link = linkConverted;
//	}
};


}


#endif /* LINKTOGAUGEFIELDCONVERT_H_ */
