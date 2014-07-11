/**
 * This is a quick and dirty solution.
 * TODO needs optimization!
 */

#ifndef GELLMANNMATRICES_H_
#define GELLMANNMATRICES_H_
#include "../../lattice/LocalLink.h"
#include "../../lattice/parameterization_types/SUNComplexFull.h"
#include "../../cuLGT1legacy/Complex.hxx"
#include "../../cudacommon/cuda_host_device.h"

using culgt::LocalLink;
using culgt::SUNComplexFull;
using culgt::Complex;


template<typename LinkType> class GellMannMatrices
{
public:
	/**
	 * We generate the matrices in SUNComplexFull Format and convert it to the wanted type in the end...
	 * @param i
	 * @return
	 */
	static CUDA_HOST_DEVICE LinkType get( int i )
	{
		typedef Complex<typename LinkType::PARAMTYPE::REALTYPE> MyComplex;
		LinkType result;
		LocalLink<SUNComplexFull<3,typename LinkType::PARAMTYPE::REALTYPE> > link;
		link.zero();
		switch( i )
		{
		case 0:
			link.identity();
			break;
		case 1:
			link.set( 1, MyComplex(1,0) );
			link.set( 3, MyComplex(1,0) );
			break;
		case 2:
			link.set( 1, MyComplex(0,-1) );
			link.set( 3, MyComplex(0,1) );
			break;
		case 3:
			link.set( 0, MyComplex(1,0) );
			link.set( 4, MyComplex(-1,0) );
			break;
		case 4:
			link.set( 2, MyComplex(1,0) );
			link.set( 6, MyComplex(1,0) );
			break;
		case 5:
			link.set( 2, MyComplex(0,-1) );
			link.set( 6, MyComplex(0,1) );
			break;
		case 6:
			link.set( 5, MyComplex(1,0) );
			link.set( 7, MyComplex(1,0) );
			break;
		case 7:
			link.set( 5, MyComplex(0,-1) );
			link.set( 7, MyComplex(0,1) );
			break;
		case 8:
			link.set( 0, MyComplex(1./sqrt(3.),0) );
			link.set( 4, MyComplex(1./sqrt(3.),0) );
			link.set( 8, MyComplex(-2./sqrt(3.),0) );
			break;
		}

		result = link;
		return result;
	}
};


#endif /* GELLMANNMATRICES_H_ */
