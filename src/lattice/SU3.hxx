/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 *
 * SU3 serves as a wrapper class for different storage techniques and implements some important operations for SU3-matrices,
 * like summation, multiplication, det(), trace(), ...
 *
 * At present, you can choose among two storage techniques:
 *  1) Link: Only information to the position of the matrix in a large array is stored. No memory is allocated!
 *  2) Matrix: Memory is allocated for storage of the SU3-matrix.
 * All operations that have two SU3-objects as input (like *, +, =, ...) have to be implement such that, they allow objects
 * with different storage technique, i.e. SU3<Link> + SU3<Matrix> should be implemented.
 * See for example the implementation of operator=().
 * All storage techniques have to implement a basic set of functions.
 *
 * TODO:
 *  - Implement all mathematical operations: like operator*
 *  - We have to introduce a compatible and clear way to switch between full 18 parameter matrices
 *    and 12 parameter (third-line-reconstruction) techniques.
 *  - At present we use the hardcodeed "Complex<Real>". This is no good style... Introduce a template parameter or something...
 */

#ifndef SU3_HXX_
#define SU3_HXX_

#include "cuda/cuda_host_device.h"
#include "datatype/datatypes.h"
#include "datatype/lattice_typedefs.h"
#include "Matrix.hxx"
#include "Complex.hxx"
#include "Link.hxx"
#include "Quaternion.hxx"
#include <math.h>

template<class Type> class SU3
{
public:
	CUDA_HOST_DEVICE SU3( Type mat );
	CUDA_HOST_DEVICE SU3();
	Type mat;

	CUDA_HOST_DEVICE inline Type& getMat();

	CUDA_HOST_DEVICE inline Complex<Real> get( lat_group_dim_t i, lat_group_dim_t j );
	CUDA_HOST_DEVICE inline Complex<Real> get(lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j);
	CUDA_HOST_DEVICE inline Quaternion<Real> getSubgroupQuaternion( lat_group_dim_t iSub, lat_group_dim_t jSub ); // TODO binding this class to class Quaternion is not good style -> make this a static function elsewhere
	CUDA_HOST_DEVICE inline Matrix<Complex<Real>, 2> getSubgroupMatrix( lat_group_dim_t iSub, lat_group_dim_t jSub );
	CUDA_HOST_DEVICE inline Matrix<Complex<Real>, 2> getSubgroupMatrixHermitian( lat_group_dim_t iSub, lat_group_dim_t jSub );
	CUDA_HOST_DEVICE inline void set( lat_group_dim_t i, lat_group_dim_t j, Complex<Real> c);
	CUDA_HOST_DEVICE inline void set(lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j, Complex<Real> c);
	CUDA_HOST_DEVICE inline SU3<Type>& operator+=( SU3<Type> ); // TODO overload for types like SU3<Link>::operator+=( SU3<Matrix> )
	CUDA_HOST_DEVICE inline SU3<Type>& operator-=( SU3<Type> ); // TODO overload for types like SU3<Link>::operator+=( SU3<Matrix> )
	CUDA_HOST_DEVICE inline SU3<Type>& operator*=( SU3<Type> );  // TODO can we overload this for arbitrary Type2 (we need some temporary matrix but we don't know of which type we want it!
																	// The problem boils down to that we can't do SU3<Link>*SU3<Link> because we don't know where to store the data...
																	// Maybe the SU3 wrapper is not as good as it looked in the first place.

	CUDA_HOST_DEVICE inline SU3<Type>& operator/=( Complex<Real> );
	CUDA_HOST_DEVICE inline SU3<Type>& operator-=( Complex<Real> );

	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& operator+=( SU3<Type2> );
	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& operator=( SU3<Type2> );
//	template<class Type2> CUDA_HOST_DEVICE inline SU3<Matrix<Complex<Real>,3> > operator*( SU3<Type2> );
	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& assignWithoutThirdLine( SU3<Type2> );
//	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& assignWithoutThirdLineFloat4ToMat( SU3<Type2> ); // TODO uebelster HACK
//	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& assignWithoutThirdLineFloat4FromMat( SU3<Type2> );// TODO uebelster HACK

	CUDA_HOST_DEVICE inline Complex<Real> det();
	CUDA_HOST_DEVICE inline Complex<Real> trace();
	CUDA_HOST_DEVICE inline void identity();
	CUDA_HOST_DEVICE inline void zero();
	CUDA_HOST_DEVICE inline void projectSU3();
	CUDA_HOST_DEVICE inline void projectSU3withoutThirdRow();
	CUDA_HOST_DEVICE inline SU3<Type>& hermitian();

	CUDA_HOST_DEVICE inline void reconstructThirdLine();

//	CUDA_HOST_DEVICE inline void leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] );
//	CUDA_HOST_DEVICE inline void rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] );

	CUDA_HOST_DEVICE inline void leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q );
	CUDA_HOST_DEVICE inline void rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q );

	CUDA_HOST_DEVICE inline void print();

	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type> operator*( SU3<Type2> b  );
	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type> operator+( SU3<Type2> b  );


	//TODO what if we compile on g++? We can't have a __device__ function!!!
#ifdef __CUDACC__
	template<class SubgroupOperationClass> __device__ static inline void perSubgroup(SubgroupOperationClass t);
#endif
};

/**
 * Creates a new SU3 object with from the given background storage object 'mat'.
 * @parameter Object that holds the information of the SU3-matrix.
 */
template<class Type> SU3<Type>::SU3( Type mat ) : mat(mat)
{
}

/**
 * TODO check if we need this.
 */
template<class Type> SU3<Type>::SU3()
{
}

template<class Type> Type& SU3<Type>::getMat()
{
	return mat;
}

/**
 * Get matrix element. Delegates the get to the underlying storage class.
 * @parameter row index i
 * @parameter col index j
 * @return matrix element (i,j)
 */
template<class Type> Complex<Real> SU3<Type>::get( lat_group_dim_t i, lat_group_dim_t j )
{
	return mat.get(i,j);
}

/**
 * Set matrix element. Delegates the set to the underlying storage class.
 * @parameter row index i
 * @parameter col index j
 * @parameter element to set
 */
template<class Type> void SU3<Type>::set( lat_group_dim_t i, lat_group_dim_t j, Complex<Real> c )
{
	return mat.set(i,j,c);
}

/**
 * Get matrix element within a 2x2-submatrix from which the SU2-subgroup is built. The indexing works as follows:
 * iSub != jSub define the 2x2-matrix A as: a_00 = M(iSub,iSub), a_01 = M(iSub,jSub), a_10 = M(jSub,iSub), a_11 = M(jSub,jSub).
 * @parameter first subgroup index
 * @parameter second subgroup index
 * @parameter row of 2x2-subgroup, i.e. i = 0 or 1
 * @parameter col of 2x2-subgroup, i.e. j = 0 or 1
 * @return matrix element in 2x2-submatrix.
 */
template<class Type> Complex<Real> SU3<Type>::get( lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j )
{
	return mat.get( (i==0)?(iSub):(jSub), (j==1)?(iSub):(jSub) );
}

/**
 * Set matrix element within a 2x2-submatrix. See corresponding get().
 * @parameter first subgroup index
 * @parameter second subgroup index
 * @parameter row of 2x2-subgroup, i.e. i = 0 or 1
 * @parameter col of 2x2-subgroup, i.e. j = 0 or 1
 * @parameter element to set.
 */
template<class Type> void SU3<Type>::set( lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j, Complex<Real> c )
{
	return mat.set((i==0)?(iSub):(jSub), (j==1)?(iSub):(jSub) ,c);
}

/*
 * Returns a 2x2-complex-subgroup in Quaternion (4 real parameter) representation. THE ELEMENT IS NOT NORMALIZED TO SU(2),
 * i.e. it is not a SU(2) element but proportional to a SU(2) element.
 * See the get()-subgroup-function for an explanation of how the subgroup is choosen.
 *
 * The function takes the 2x2-matrix and mixes the corresponding elements, see for example Gattringer & Lang (2010), p. 81.
 *
 * TODO imply non-normalization by function name.
 *
 * @parameter first subgroup index
 * @parameter second subgroup index
 */
template<class Type> Quaternion<Real> SU3<Type>::getSubgroupQuaternion( lat_group_dim_t iSub, lat_group_dim_t jSub )
{
	Quaternion<Real> q;
	Complex<Real> temp;
	temp = mat.get(iSub,iSub);
	q[0] = temp.x;
	q[3] = temp.y;
	temp = mat.get(jSub,jSub);
	q[0] += temp.x;
	q[3] -= temp.y;
	temp = mat.get(iSub,jSub);
	q[2] = temp.x;
	q[1] = temp.y;
	temp = mat.get(jSub,iSub);
	q[2] -= temp.x;
	q[1] += temp.y;


//	q[0] *= .5;
//	q[1] *= .5;
//	q[2] *= .5;
//	q[3] *= .5;

	return q;
}

template<class Type> Matrix<Complex<Real>, 2 > SU3<Type>::getSubgroupMatrix( lat_group_dim_t iSub, lat_group_dim_t jSub )
{
	Matrix<Complex<Real>, 2 > subMatrix;

	subMatrix.set( 0, 0, mat.get(iSub, iSub) );
	subMatrix.set( 0, 1, mat.get(iSub, jSub) );
	subMatrix.set( 1, 0, mat.get(jSub, iSub) );
	subMatrix.set( 1, 1, mat.get(jSub, jSub) );

	return subMatrix;
}

template<class Type> Matrix<Complex<Real>, 2 > SU3<Type>::getSubgroupMatrixHermitian( lat_group_dim_t iSub, lat_group_dim_t jSub )
{
	Matrix<Complex<Real>, 2 > subMatrix;

	subMatrix.set( 0, 0, mat.get(iSub, iSub).conj() );
	subMatrix.set( 0, 1, mat.get(jSub, iSub).conj() );
	subMatrix.set( 1, 0, mat.get(iSub, jSub).conj() );
	subMatrix.set( 1, 1, mat.get(jSub, jSub).conj() );

	return subMatrix;
}

/**
 * Delegate operation to underlying storage class.
 */
template<class Type> SU3<Type>& SU3<Type>::operator+=( SU3<Type> c )
{
	mat += c.mat;
	return *this;
}

/**
 *
 */
template<class Type> template<class Type2> SU3<Type>& SU3<Type>::operator+=( SU3<Type2> c )
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			mat.set(i,j,mat.get(i,j) + c.mat.get(i,j));
		}
	return *this;
}

/**
 * Delegate operation to underlying storage class.
 */
template<class Type> SU3<Type>& SU3<Type>::operator-=( SU3<Type> c )
{
	mat -= c.mat;
	return *this;
}

/**
 * Delegate operation to underlying storage class.
 */
template<class Type> SU3<Type>& SU3<Type>::operator*=( SU3<Type> c )
{
	mat *= c.mat;
	return *this;
}

template<class Type> SU3<Type>& SU3<Type>::operator/=( Complex<Real> c )
{
	mat /= c;
	return *this;
}

template<class Type> SU3<Type>& SU3<Type>::operator-=( Complex<Real> c )
{
	mat -= c;
	return *this;
}

/**
 * Assignment operator allowing different storage class at right and left side of the assignement.
 */
template<class Type> template<class Type2> SU3<Type>& SU3<Type>::operator=( SU3<Type2> c )
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			mat.set(i,j,c.mat.get(i,j));
		}
	return *this;
}

///**
// * TODO comments!
// * We return always with storage type Matrix. Think about how this can be done in a nice way.
// */
//template<class Type2> CUDA_HOST_DEVICE inline SU3<Matrix<Complex<Real>,3> > operator*( SU3<Type2> )
//{
//
//}

/**
 * Assignement like operator=, but does not copy the third line.
 */
template<class Type> template<class Type2> SU3<Type>& SU3<Type>::assignWithoutThirdLine( SU3<Type2> c )
{
	for( lat_group_dim_t i = 0; i < 2; i++ )
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			mat.set(i,j,c.mat.get(i,j));
		}
	return *this;
}


/**
 * TODO for test: assume type2==link, type==matrix
 */
//template<class Type> template<class Type2> SU3<Type>& SU3<Type>::assignWithoutThirdLineFloat4ToMat( SU3<Type2> c )
//{
//	float4 f0,f1,f2;
//	f0 = c.mat.getFloat4( 0, 0 );
//	f1 = c.mat.getFloat4( 0, 2 );
//	f2 = c.mat.getFloat4( 1, 1 );
//
//	Complex<Real> a( f0.x, f0.y );
//	mat.set( 0, 0, a);
//
//	a = Complex<Real>( f0.z, f0.w );
//	mat.set( 0, 1, a);
//
//	a = Complex<Real>( f1.x, f1.y );
//		mat.set( 0, 2, a);
//
//	a = Complex<Real>( f1.z, f1.w );
//		mat.set( 1, 0, a);
//
//	a = Complex<Real>( f2.x, f2.y );
//	mat.set( 1, 1, a);
//
//	a = Complex<Real>( f2.z, f2.w );
//	mat.set( 1, 2, a);
//
//	return *this;
//}
//
//template<class Type> template<class Type2> SU3<Type>& SU3<Type>::assignWithoutThirdLineFloat4FromMat( SU3<Type2> c )
//{
//	float4 f0,f1,f2;
//
//	f0.x = c.mat.get(0,0).x;
//	f0.y = c.mat.get(0,0).y;
//	f0.z = c.mat.get(0,1).x;
//	f0.w = c.mat.get(0,1).y;
//
//	f1.x = c.mat.get(0,2).x;
//	f1.y = c.mat.get(0,2).y;
//	f1.z = c.mat.get(1,0).x;
//	f1.w = c.mat.get(1,0).y;
//
//	f2.x = c.mat.get(1,1).x;
//	f2.y = c.mat.get(1,1).y;
//	f2.z = c.mat.get(1,2).x;
//	f2.w = c.mat.get(1,2).y;
//
//
//	mat.setFloat4( 0, 0, f0 );
//	mat.setFloat4( 0, 2, f1 );
//	mat.setFloat4( 1, 1, f2 );
//	return *this;
//}

/**
 * Multiplication of SU3 matrices.
 */
template<class Type> template<class Type2> SU3<Type> SU3<Type>::operator*( SU3<Type2> b )
{
	SU3<Type> c;
	c.zero();
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			Complex<Real> temp(0,0);
			for( lat_group_dim_t k = 0; k < 3; k++ )
			{
				temp += get(i,k)*b.get(k,j);
			}
			c.set(i,j,temp);
		}
	}
	return c;
}

/**
 * Summation of SU3 matrices.
 */
template<class Type> template<class Type2> SU3<Type> SU3<Type>::operator+( SU3<Type2> b )
{
	SU3<Type> c(*this);
	return c+=b;
}

/**
 * Determinant.
 * TODO Check for faster implementations...
 */
template<class Type> Complex<Real> SU3<Type>::det()
{
	Complex<Real> c( 0, 0 );
	Complex<Real> temp = get(0,0);
	temp *= get(1,1);
	temp *= get(2,2);
	c += temp;

	temp = get(0,1);
	temp *= get(1,2);
	temp *= get(2,0);
	c += temp;

	temp = get(0,2);
	temp *= get(1,0);
	temp *= get(2,1);
	c += temp;



	temp = get(2,0);
	temp *= get(1,1);
	temp *= get(0,2);
	c -= temp;

	temp = get(1,0);
	temp *= get(0,1);
	temp *= get(2,2);
	c -= temp;

	temp = get(0,0);
	temp *= get(2,1);
	temp *= get(1,2);
	c -= temp;

	return c;
}

/**
 * Delegate trace() to underlying storage class.
 */
template<class Type> Complex<Real> SU3<Type>::trace()
{
	return mat.trace();
}

/**
 * Set matrix to identity.
 */
template<class Type> void SU3<Type>::identity()
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for(lat_group_dim_t j = 0; j < 3; j++ )
		{
			if( i == j )
			{
				set( i,j, Complex<Real>( 1, 0 ) );
			}
			else
			{
				set( i,j, Complex<Real>( 0, 0 ) );
			}
		}
	}
}

/**
 * Set matrix to its hermitian. TODO this is a straight forward implementation, do it in a fast way!!!
 */
template<class Type> SU3<Type>& SU3<Type>::hermitian()
{
	Matrix<Complex<Real>,3> temp;
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			temp.set(j,i,get(i,j).conj());
		}
	}
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			set(i,j,temp.get(i,j));
		}
	}
	return *this;
}

/**
 * Set matrix to zero.
 */
template<class Type> void SU3<Type>::zero()
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for(lat_group_dim_t j = 0; j < 3; j++ )
		{
			set( i,j, Complex<Real>(0, 0 ) );
		}
	}
}

/**
 * Project back to SU3. Needed because of numerical inaccuracy.
 *
 * Technique:
 * - Normalize fist row
 * - Orthogonalize second row with respect to the first row.
 * - Normalize second row
 * - Reconstruct the third row from the first two rows (cross-product).
 */
template<class Type> void SU3<Type>::projectSU3()
{
	Real abs_u = 0, abs_v = 0;
	Complex<Real> sp(0.,0.);

	// normalize first row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		abs_u += get( 0, i ).abs_squared();
	}

	abs_u = sqrt(abs_u);

	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 0, i, get(0,i)/abs_u );
	}

	// orthogonalize second row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		sp += get( 1,i ) * get( 0, i ).conj();
	}
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 1, i, get(1,i) - get( 0, i)*sp );
	}

	// normalize second row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		abs_v += get( 1, i ).abs_squared();
	}
	abs_v = sqrt(abs_v);
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 1, i, get(1,i)/abs_v );
	}

	// reconstruct third row
	reconstructThirdLine();

//	set( 2, 0, get(0,1).conj() * get(1,2).conj() - get(0,2).conj() * get(1,1).conj() );
//
//	set( 2, 1, get(0,2).conj() * get(1,0).conj() - get(0,0).conj() * get(1,2).conj() );
//
//	set( 2, 2, get(0,0).conj() * get(1,1).conj() - get(0,1).conj() * get(1,0).conj() );
}

/**
 * Project back to SU3. Needed because of numerical inaccuracy.
 * Does not touch third row. TODO: We have to find a way to do this implicit.
 *
 * Technique:
 * - Normalize fist row
 * - Orthogonalize second row with respect to the first row.
 * - Normalize second row
 * - Reconstruct the third row from the first two rows (cross-product).
 */
template<class Type> void SU3<Type>::projectSU3withoutThirdRow()
{
	Real abs_u = 0, abs_v = 0;
	Complex<Real> sp(0.,0.);

	// normalize first row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		abs_u += get( 0, i ).abs_squared();
	}

	abs_u = sqrt(abs_u);

	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 0, i, get(0,i)/abs_u );
	}

	// orthogonalize second row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		sp += get( 1,i ) * get( 0, i ).conj();
	}
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 1, i, get(1,i) - get( 0, i)*sp );
	}

	// normalize second row
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		abs_v += get( 1, i ).abs_squared();
	}
	abs_v = sqrt(abs_v);
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		set( 1, i, get(1,i)/abs_v );
	}
}

/**
 * Reconstructs the third line of the 3x3 matrix from the first two rows. This is needed when we load only the first two rows
 * from memory to save memory bandwidth in CUDA applications.
 */
template<class Type> void SU3<Type>::reconstructThirdLine()
{
//	Complex<Real> a00 = get(0,0);
//	Complex<Real> a01 = get(0,1);
//	Complex<Real> a02 = get(0,2);
//	Complex<Real> a10 = get(1,0);
//	Complex<Real> a11 = get(1,1);
//	Complex<Real> a12 = get(1,2);
//
//	set( 2, 0, Complex<Real>( a01.x*a12.x-a01.y*a12.y - a11.x*a02.x+a11.y*a02.y, -a01.x*a12.y-a01.y*a12.x + a11.x*a02.y+a11.y*a02.x ) );
//	set( 2, 1, Complex<Real>( a02.x*a10.x-a02.y*a10.y - a12.x*a00.x+a12.y*a00.y, -a02.x*a10.y-a02.y*a10.x + a12.x*a00.y+a12.y*a00.x ) );
//	set( 2, 2, Complex<Real>( a00.x*a11.x-a00.y*a11.y - a01.x*a10.x+a01.y*a10.y, -a00.x*a11.y-a00.y*a11.x + a01.x*a10.y+a01.y*a10.x ) );




	Complex<Real> a = get(0,1).conj() * get(1,2).conj() - get(0,2).conj() * get(1,1).conj();
	Complex<Real> b = get(0,2).conj() * get(1,0).conj() - get(0,0).conj() * get(1,2).conj();
	Complex<Real> c = get(0,0).conj() * get(1,1).conj() - get(0,1).conj() * get(1,0).conj();

	Real norm = a.abs_squared() + b.abs_squared() + c.abs_squared();
	norm = 1./sqrt( norm );

	set( 2, 0, a*norm );
	set( 2, 1, b*norm );
	set( 2, 2, c*norm );


//	Complex<Real> a = get(0,1).conj() * get(1,2).conj() - get(0,2).conj() * get(1,1).conj();
//	Complex<Real> b = get(0,2).conj() * get(1,0).conj() - get(0,0).conj() * get(1,2).conj();
//	Complex<Real> c = get(0,0).conj() * get(1,1).conj() - get(0,1).conj() * get(1,0).conj();
//
//	set( 2, 0, a);
//
//	set( 2, 1, b);
//
//	set( 2, 2, c );
}

/**
 * Muliplication of the SU3 matrix by a SU2-subgroup element in Quaternion representation from the right.
 * TODO Be a bit more precise in the commentary...
 */
template<class Type> void SU3<Type>::rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q )
{
	for( lat_group_dim_t k = 0; k < 3; k++ )
	{
//		ki = k*3+iSub;
		Complex<Real> KI = q->get( 0, 0 ) * get(k,i);
		KI += q->get( 1, 0 ) * get(k,j);

		Complex<Real> KJ = q->get( 0, 1 ) * get(k,i);
		KJ += q->get(1,1) * get(k,j);

		set( k, i , KI);
		set( k, j , KJ);
	}
}

/**
 * Muliplication of the SU3 matrix by a SU2-subgroup element in Quaternion representation from the left.
 */
template<class Type> void SU3<Type>::leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q )
{
	for( lat_group_dim_t k = 0; k < 3; k++ )
	{
		Complex<Real> IK = q->get( 0, 0 ) * get(i,k);
		IK += q->get( 0, 1 ) * get(j,k);

		Complex<Real> JK = q->get( 1, 0 ) * get(i,k);
		JK += q->get(1,1) * get(j,k);

		set( i, k , IK );
		set( j, k,  JK );
	}
}

template<class Type> void SU3<Type>::print()
{
	printf( "[%f+i*%f\t %f+i*%f\t %f+i*%f\n", get(0,0).x, get(0,0).y, get(0,1).x, get(0,1).y, get(0,2).x, get(0,2).y );
	printf( "%f+i*%f\t %f+i*%f\t %f+i*%f\n", get(1,0).x, get(1,0).y, get(1,1).x, get(1,1).y, get(1,2).x, get(1,2).y );
	printf( "%f+i*%f\t %f+i*%f\t %f+i*%f]\n", get(2,0).x, get(2,0).y, get(2,1).x, get(2,1).y, get(2,2).x, get(2,2).y );
}


//template<class Type> void SU3<Type>::rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] )
//{
//	for( lat_group_dim_t k = 0; k < 3; k++ )
//	{
////		ki = k*3+iSub;
//		Complex<Real> KI = Complex<Real>(q[0],q[3]) * get(k,i);
//		KI += Complex<Real>(-q[2],q[1]) * get(k,j);
//
//		Complex<Real> KJ = Complex<Real>(q[2],q[1]) * get(k,i);
//		KJ += Complex<Real>(q[0],-q[3]) * get(k,j);
//
//		set( k, i , KI);
//		set( k, j , KJ);
//	}
//}
//
//template<class Type> void SU3<Type>::leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] )
//{
//	for( lat_group_dim_t k = 0; k < 3; k++ )
//	{
//		Complex<Real> IK = Complex<Real>(q[0],q[3]) * get(i,k);
//		IK += Complex<Real>(q[2],q[1])* get(j,k);
//
//		Complex<Real> JK = Complex<Real>(-q[2], q[1]) * get(i,k);
//		JK += Complex<Real>(q[0],-q[3]) * get(j,k);
//
//		set( i, k , IK );
//		set( j, k,  JK );
//	}
//}



#ifdef __CUDACC__
/**
 * Performs an operation defined in the SubgroupOperationClass by calling its subgroup() function for 3 SU3 subgroups.
 * Check an example application, like the Coulomb-gaugefixing routine.
 */
template<class Type> template<class SubgroupOperationClass> __device__ void SU3<Type>::perSubgroup( SubgroupOperationClass t )
{
//	t.subgroup(0,1);
//	t.subgroup(0,2);
//	t.subgroup(1,2);

	t.subgroup(0,2);
	t.subgroup(1,2);
	t.subgroup(0,1);
}
#endif


#endif /* SU3_HXX_ */
