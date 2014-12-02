/**
 * matrix_vector.h
 *
 *  Created on: Jul 7, 2014
 *      Author: vogt
 */

#ifndef MATRIX_VECTOR_H_
#define MATRIX_VECTOR_H_

#include "Vector.h"
#include "Matrix.h"

namespace culgt
{


template<typename T, int Size> CUDA_HOST_DEVICE  Vector<T,Size> operator*( Matrix<T,Size>& lhs, Vector<T,Size>& rhs )
{
	Vector<T,Size> result;
	for( int i = 0; i < Size; i++ )
	{
		result(i) = 0.;
		for( int j = 0; j < Size; j++ )
		{
			result(i) += lhs(i,j)*rhs(j);
		}
	}
	return result;
}

}
#endif /* MATRIX_VECTOR_H_ */
