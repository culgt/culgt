/**
 * Storage class for a lattice link stored in a big array that keeps the complete configuration. The array has to allocated elsewhere.
 * Only the access to the matrix-elements is computed and a pointer to the array is kept.
 *
 * TODO:
 *  - Why is T_Ndim a template parameter? All information is encoded in TheSite? Remove this.
 *  - Do we want a minimal set of functions here, like only getter and setter OR do we wannt to define other operations here like +, *, ...
 *    The question is: Is it possible to write mor efficient functions in the storage classes rather than in the frontend-classes like "SU3"?
 *  - Do not use typedef'ed "complex".
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-04-16
 */

#ifndef LINK_HXX_
#define LINK_HXX_

#include <assert.h>
#include <iostream>
#include "datatype/datatypes.h"
#include "Complex.hxx"

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> class Link
{
public:
	CUDA_HOST_DEVICE inline Link( Real* data, TheSite site, int mu );
	CUDA_HOST_DEVICE inline virtual ~Link();
	CUDA_HOST_DEVICE inline Complex<Real> get(int i, int j);
//	CUDA_HOST_DEVICE inline float4 getFloat4(int i, int j);
	CUDA_HOST_DEVICE inline void set(int i, int j, Complex<Real> c);
//	CUDA_HOST_DEVICE inline void setFloat4(int i, int j, float4 f);
	CUDA_HOST_DEVICE inline Complex<Real> trace();
	CUDA_HOST_DEVICE inline Link<Pattern, TheSite, T_Ndim, T_Nc>& operator+=( Link<Pattern, TheSite, T_Ndim, T_Nc> );

	CUDA_HOST_DEVICE inline TheSite& getSite();
	CUDA_HOST_DEVICE inline void setMu( int mu );
	CUDA_HOST_DEVICE inline void setPointer( Real* pointer );

private:
	Real* data; // pointer to the link array
	TheSite site; // current lattice site
	int mu; // direction of the link
};

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Link<Pattern, TheSite, T_Ndim, T_Nc>::Link( Real* data, TheSite site, int mu ) : data(data), site( site ), mu(mu)
{
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Link<Pattern, TheSite, T_Ndim, T_Nc>::~Link()
{
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> TheSite& Link<Pattern, TheSite, T_Ndim, T_Nc>::getSite()
{
	return site;
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> void Link<Pattern, TheSite, T_Ndim, T_Nc>::setMu( int mu )
{
	this->mu = mu;
}

template<class Pattern, class TheSite, int T_Ndim, int T_Nc> void Link<Pattern, TheSite, T_Ndim, T_Nc>::setPointer( Real* pointer )
{
	this->data = pointer;
}

//template<class Pattern, class TheSite, int T_Ndim, int T_Nc> float4 Link<Pattern, TheSite, T_Ndim, T_Nc>::getFloat4( int i, int j )
//{
//	float4 *fdata = (float4*)data;
//	return fdata[Pattern::getIndex( site, mu, i, j, 0 )/4];
//}
//
//template<class Pattern, class TheSite, int T_Ndim, int T_Nc> void Link<Pattern, TheSite, T_Ndim, T_Nc>::setFloat4( int i, int j, float4 f )
//{
//	float4 *fdata = (float4*)data;
//	fdata[Pattern::getIndex( site, mu, i, j, 0 )/4] = f;
//
//}

/**
 * Returns the matrix element (i,j).
 * @parameter row index i
 * @parameter col index j
 * @return element (i,j)
 */
template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Complex<Real> Link<Pattern, TheSite, T_Ndim, T_Nc>::get( int i, int j )
{
	return Complex<Real>( data[Pattern::getIndex( site, mu, i, j, 0 )], data[Pattern::getIndex( site, mu, i, j, 1 )] );

}

/**
 * Sets the matrix element (i,j).
 * @parameter row index i
 * @parameter col index j
 * @parameter element to set
 */
template<class Pattern, class TheSite, int T_Ndim, int T_Nc> void Link<Pattern, TheSite, T_Ndim, T_Nc>::set( int i, int j, Complex<Real> c )
{
	data[Pattern::getIndex( site, mu, i, j, 0 )] = c.x;
	data[Pattern::getIndex( site, mu, i, j, 1 )] = c.y;
}

/**
 * Trace.
 * TODO do it here or in frontend class SU3? Maybe we want to use Link without frontend class?
 */
template<class Pattern, class TheSite, int T_Ndim, int T_Nc> Complex<Real> Link<Pattern, TheSite, T_Ndim, T_Nc>::trace()
{
	Complex<Real> c;
	for( int i = 0; i < T_Nc; i++ )
	{
		c += get(i,i);
	}
	return c;
}

/**
 * Add and assign...
 */
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
