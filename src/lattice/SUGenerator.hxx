/*
 * SUGenerator.hxx
 *
 *  Created on: Jul 3, 2012
 *      Author: vogt
 *
 * TODO think about storing products T^a*T^b and {T^a,T^b} and things like this...
 * TODO I did not test this class in CUDA code. Maybe we don't want to have precalculated products in CUDA
 * 		OR we want to have them in constant memory!
 */

#ifndef SUGENERATOR_HXX_
#define SUGENERATOR_HXX_

#include "../util/datatype/lattice_typedefs.h"
#include "../util/datatype/datatypes.h"
#include "Matrix.hxx"
#include "assert.h"
#include <math.h>

template<int Nc> class SUGenerator
{
public:
	static void init();
	static void initProduct();
	static void initAntiCommutator();
	static Matrix<complex,Nc>& get( int i );
	static Matrix<complex,Nc>& getProduct( int a, int b );
	static Matrix<complex,Nc>& getAntiCommutator( int a, int b );
private:
	static bool initialized;
	static bool initializedProduct;
	static bool initializedAntiCommutator;
	static Matrix<complex,Nc> generator[Nc*Nc];
	static Matrix<complex,Nc> product[(Nc*Nc-1)*(Nc*Nc-1)];
	static Matrix<complex,Nc> antiCommutator[(Nc*Nc-1)*(Nc*Nc-1)];
};

template<int Nc> bool SUGenerator<Nc>::initialized = false;
template<int Nc> bool SUGenerator<Nc>::initializedProduct = false;
template<int Nc> bool SUGenerator<Nc>::initializedAntiCommutator = false;
template<int Nc> Matrix<complex,Nc> SUGenerator<Nc>::generator[Nc*Nc];
template<int Nc> Matrix<complex,Nc> SUGenerator<Nc>::product[(Nc*Nc-1)*(Nc*Nc-1)];
template<int Nc> Matrix<complex,Nc> SUGenerator<Nc>::antiCommutator[(Nc*Nc-1)*(Nc*Nc-1)];

template<> void SUGenerator<2>::init()
{
	initialized = true;

	generator[0] = complex(1,0);

	generator[1] = complex(0,0);
	generator[1].set( 0, 1, complex( 0, 1 ) );
	generator[1].set( 1, 0, complex( 0, 1 ) );

	generator[2] = complex(0,0);
	generator[2].set( 0, 1, complex( -1, 0 ) );
	generator[2].set( 1, 0, complex( 1, 0 ) );

	generator[3] = complex(0,0);
	generator[3].set( 0, 0, complex( 0, 1 ) );
	generator[3].set( 1, 1, complex( 0, -1 ) );
}

template<> void SUGenerator<3>::init()
{
	initialized = true;

	generator[0] = complex(1,0);

	generator[1] = complex(0,0);
	generator[1].set( 0, 1, complex( 1, 0 ) );
	generator[1].set( 1, 0, complex( 1, 0 ) );

	generator[2] = complex(0,0);
	generator[2].set( 0, 1, complex( 0, -1 ) );
	generator[2].set( 1, 0, complex( 0, 1 ) );

	generator[3] = complex(0,0);
	generator[3].set( 0, 0, complex( 1, 0 ) );
	generator[3].set( 1, 1, complex( -1, 0 ) );

	generator[4] = complex(0,0);
	generator[4].set( 0, 2, complex( 1, 0 ) );
	generator[4].set( 2, 0, complex( 1, 0 ) );

	generator[5] = complex(0,0);
	generator[5].set( 0, 2, complex( 0, -1 ) );
	generator[5].set( 2, 0, complex( 0, 1 ) );

	generator[6] = complex(0,0);
	generator[6].set( 1, 2, complex( 1, 0 ) );
	generator[6].set( 2, 1, complex( 1, 0 ) );

	generator[7] = complex(0,0);
	generator[7].set( 1, 2, complex( 0, -1 ) );
	generator[7].set( 2, 1, complex( 0, 1 ) );

	generator[8] = complex(0,0);
	generator[8].set( 0, 0, complex( 1, 0 ) );
	generator[8].set( 1, 1, complex( 1, 0 ) );
	generator[8].set( 2, 2, complex( -2, 0 ) );
	generator[8] *= 1./pow(3.,.5);


}

template<int Nc> void SUGenerator<Nc>::init()
{
	// TODO implement for arbitary Nc, see for example fermiqcd_su_generator
	assert(false);
}

template<int Nc> void SUGenerator<Nc>::initProduct()
{
	if( initialized ) // The generators have to be initialized
	{
		initializedProduct = true;
		for( int a = 0; a < Nc*Nc-1; a++ )
			for( int b = 0; b < Nc*Nc-1; b++ )
			{
				product[a*(Nc*Nc-1)+b] = generator[a+1];
				product[a*(Nc*Nc-1)+b] *= generator[b+1];
			}
	}
}

template<int Nc> void SUGenerator<Nc>::initAntiCommutator()
{
	if( initializedProduct ) // The generators have to be initialized
	{
		initializedAntiCommutator = true;
		for( int a = 0; a < Nc*Nc-1; a++ )
			for( int b = 0; b < Nc*Nc-1; b++ )
			{
				antiCommutator[a*(Nc*Nc-1)+b] = product[a*(Nc*Nc-1)+b];
				antiCommutator[a*(Nc*Nc-1)+b] += product[b*(Nc*Nc-1)+a];
			}
	}
}

template<int Nc> Matrix<complex,Nc>& SUGenerator<Nc>::get( int i )
{
	if( initialized ) return generator[i];
	else assert( false ); // TODO put some nice error handling here
}

template<int Nc> Matrix<complex,Nc>& SUGenerator<Nc>::getProduct( int a, int b )
{
	if( initializedProduct ) return product[(a-1)*(Nc*Nc-1)+(b-1)];
	else assert( false ); // TODO put some nice error handling here
}

template<int Nc> Matrix<complex,Nc>& SUGenerator<Nc>::getAntiCommutator( int a, int b )
{
	if( initializedAntiCommutator ) return antiCommutator[(a-1)*(Nc*Nc-1)+(b-1)];
	else assert( false ); // TODO put some nice error handling here
}



#endif /* SUGENERATOR_HXX_ */
