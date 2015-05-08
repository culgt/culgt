/**
 * This is a quick and dirty solution.
 * TODO needs optimization!
 */

#ifndef PAULIMATRICES_H_
#define PAULIMATRICES_H_
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SUNComplexFull.h"
#include "Complex.h"
#include "cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename LinkType> class PauliMatrices
{
public:
	/**
	 * We generate the matrices in SUNRealFull Format and convert it to the wanted type in the end...
	 * @param i
	 * @return
	 */
	static CUDA_HOST_DEVICE LinkType get( int i )
	{
		LinkType result;
		LocalLink<SUNRealFull<2,typename LinkType::PARAMTYPE::REALTYPE> > link;
		link.zero();
		switch( i )
		{
		case 0:
			link.identity();
			break;
		case 1:
			link.set( 2, 1. );
			link.set( 4, 1. );
			break;
		case 2:
			link.set( 3, -1. );
			link.set( 5, 1. );
			break;
		case 3:
			link.set( 0, 1 );
			link.set( 6, -1. );
			break;
		}

		result = link;
		return result;
	}
};

}
#endif
