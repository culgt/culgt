/*
 * datatypes.h
 *
 *  Created on: Apr 16, 2012
 *      Author: vogt
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#ifdef DOUBLEPRECISION
	typedef double Real;
	#define MPI_Real MPI_DOUBLE
#else
	typedef float Real;
	#define MPI_Real MPI_FLOAT
#endif
	
#endif /* DATATYPES_H_ */
