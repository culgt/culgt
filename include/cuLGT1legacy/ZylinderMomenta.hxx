/*
 * TODO mir scheint das hat ein Vollidiot programmiert
 *
 * Assume N1=N2=N3[=N...] = N
 * ZylinderMomenta.hxx
 *
 * provides zylindercut momenta (i+n,i,i,...) with n in [-cut,...,cut]
 * and starting with (1,1,0)
 *
 *  Created on: Jul 9, 2012
 *      Author: vogt
 *
 */

#ifndef ZYLINDERMOMENTA_HXX_
#define ZYLINDERMOMENTA_HXX_

#include <iostream>
#include <vector>
#include <boost/array.hpp>
#include <math.h>


template<int Ndim, typename T> class ZylinderMomenta
{
public:
	ZylinderMomenta( int Ns, int cut, bool onePerOffDiag = false );
	int size;
	std::vector<boost::array<int,Ndim> > momenta;
	T abs( int i );
	T absLattice( int i );
	void addMomentum( int (&newmom)[Ndim] );
	std::vector<boost::array<int,Ndim> >& operator[]( int i );
private:
	void init();
	int cut;
	int Ns;
	bool onePerOffDiag; // select only (d+i,d,d) not (d,d+i,d) and (d,d,d+i)
};

template<int Ndim, typename T> ZylinderMomenta<Ndim,T>::ZylinderMomenta( int Ns, int cut, bool onePerOffDiag ) : cut(cut), Ns(Ns), momenta((2*cut+1)*Ndim*Ns/2), onePerOffDiag(onePerOffDiag)
{
	init();
}

template<int Ndim, typename T> void ZylinderMomenta<Ndim,T>::init()
{
	int counter = 0;
	for( int i = 1; i <= Ns/2; i++ )
	{
		for( int c = -cut; c <= cut; c++ )
		{
			if( (i + c <= Ns/2) && (i + c ) >= 0 )
			{
				if( onePerOffDiag )
				{
					for( int j = 0; j < Ndim; j++ )
					{
						if( j == 0 )
						{
							momenta[counter][j] = i+c;
						}
						else
						{
							momenta[counter][j] = i;
						}
					}
					counter++;
				}
				else
				{
					for( int d = 0; d < Ndim; d++ )
					{
						for( int j = 0; j < Ndim; j++ )
						{
							if( j == d )
							{
								momenta[counter][j] = i+c;
							}
							else
							{
								momenta[counter][j] = i;
							}
						}
						counter++;
						if( c == 0 ) break;
					}
				}
			}
		}
	}
	size = counter;
}

template<int Ndim, typename T> T ZylinderMomenta<Ndim,T>::abs( int i )
{
	T temp = 0;
	for( int j = 0; j < Ndim; j++ )
	{
		temp += (T)momenta[i][j]*(T)momenta[i][j];
	}
	return sqrt( temp );
}

template<int Ndim, typename T> T ZylinderMomenta<Ndim,T>::absLattice( int i )
{
	T temp = 0;
	for( int j = 0; j < Ndim; j++ )
	{
		T temp2 = 2.*sin( (T)momenta[i][j] * acos(-1.) / (T)Ns );
		temp += temp2*temp2;
	}
	return sqrt( temp );
}

template<int Ndim, typename T> void ZylinderMomenta<Ndim,T>::addMomentum( int (&newmom)[Ndim] )
{
	for( int i = 0; i < Ndim; i++ )
		momenta[size][i] = newmom[i];
	size++;
}

#endif /* ZYLINDERMOMENTA_HXX_ */
