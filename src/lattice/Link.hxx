/*
 * Link.hxx
 *
 * TODO why is T_Ndim a template parameter? All information is encoded in TheSite???
 *
 *  Created on: Apr 17, 2012
 *      Author: vogt
 */

#ifndef LINK_HXX_
#define LINK_HXX_

#include "../util/datatype/datatypes.h"
#include <assert.h>
#include <iostream>

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> class Link
{
public:
	CUDA_HOST_DEVICE inline Link( Real* data, TheSite site, int mu );
	CUDA_HOST_DEVICE inline virtual ~Link();
	CUDA_HOST_DEVICE inline complex get(int i, int j);
//	CUDA_HOST_DEVICE inline complex get(int iSub, int jSub, int i, int j);
	CUDA_HOST_DEVICE inline void set(int i, int j, complex c);
//	CUDA_HOST_DEVICE inline void set(int iSub, int jSub, int i, int j, complex c);
	CUDA_HOST_DEVICE inline complex trace();
//	CUDA_HOST_DEVICE inline complex det();
	CUDA_HOST_DEVICE inline Link<Pattern, TheSite, T_Ndim, T_Nc>& operator+=( Link<Pattern, TheSite, T_Ndim, T_Nc> );
	Real* data;
private:
	TheSite site;
	int mu;
};

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Link<Pattern, TheSite, T_Ndim, T_Nc>::Link( Real* data, TheSite site, int mu ) : data(data), site( site ), mu(mu)
{

//	this->site = site;
//	this->data = data;
//	this->mu = mu;

//	std::cout << "mu " << mu << std::endl;
//	for( int i = 0; i < site.Ndim; i++ )
//	{
//		std::cout << "size " << i << ": "<< site.size[i]<<  std::endl;
//	}
//	for( int i = 0; i < site.Ndim; i++ )
//	{
//		std::cout << "site " << i << ": "<< site[i]<<  std::endl;
//	}
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Link<Pattern, TheSite, T_Ndim, T_Nc>::~Link()
{
}


template<class Pattern, class TheSite, int T_Ndim, int T_Nc> complex Link<Pattern, TheSite, T_Ndim, T_Nc>::get( int i, int j )
{
	return complex( data[Pattern::getIndex( site, mu, i, j, 0 )], data[Pattern::getIndex( site, mu, i, j, 1 )] );

}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> void Link<Pattern, TheSite, T_Ndim, T_Nc>::set( int i, int j, complex c )
{
	data[Pattern::getIndex( site, mu, i, j, 0 )] = c.x;
	data[Pattern::getIndex( site, mu, i, j, 1 )] = c.y;
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> complex Link<Pattern, TheSite, T_Ndim, T_Nc>::trace()
{
	complex c;
	for( int i = 0; i < T_Nc; i++ )
	{
		c += get(i,i);
	}
	return c;
}

//template<class Pattern, class TheSite, int T_Ndim, int T_Nc> complex Link<Pattern, TheSite, T_Ndim, T_Nc>::det()
//{
//	assert( T_Nc == 3 ); // TODO do it properly
//	complex c( 0, 0 );
//	complex temp( 1, 0 );
//	temp *= get(0,0);
//	temp *= get(1,1);
//	temp *= get(2,2);
//	c += temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,1);
//	temp *= get(1,2);
//	temp *= get(2,0);
//	c += temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,2);
//	temp *= get(1,0);
//	temp *= get(2,1);
//	c += temp;
//
//
//
//	temp = complex( 1, 0 );
//	temp *= get(2,0);
//	temp *= get(1,1);
//	temp *= get(0,2);
//	c -= temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(1,0);
//	temp *= get(0,1);
//	temp *= get(2,2);
//	c -= temp;
//
//	temp = complex( 1, 0 );
//	temp *= get(0,0);
//	temp *= get(2,1);
//	temp *= get(1,2);
//	c -= temp;
//
//	return c;
//}


template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Link<Pattern, TheSite, T_Ndim, T_Nc>& Link<Pattern, TheSite, T_Ndim, T_Nc>::operator+=( Link<Pattern, TheSite, T_Ndim, T_Nc> a )
{
	for(int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			data[Pattern::getIndex( site, mu, i, j, 0 )] += a.get(i,j).x;
			data[Pattern::getIndex( site, mu, i, j, 1 )] += a.get(i,j).y;
		}
	}
	return *this;
}



#endif /* LINK_HXX_ */
