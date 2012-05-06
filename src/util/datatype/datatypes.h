/*
 * datatypes.h
 *
 *  Created on: Apr 16, 2012
 *      Author: vogt
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include "Complex.hxx"

#ifdef DOUBLEPRECISION
	typedef double Real;
#else
	typedef float Real;
#endif

typedef Complex<Real> complex;

#endif /* DATATYPES_H_ */
