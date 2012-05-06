/*
 * A lot TODO here: overload all operators (like *...)
 * At some point we have to introduce reconstruction patterns...
 */

#ifndef SU3_HXX_
#define SU3_HXX_

#include "../util/datatype/datatypes.h"
#include "Matrix.hxx"
#include "Link.hxx"
#include "Quaternion.hxx"

/**
 * "Type" can be Link (which acts on a global arrary of SU3 matrices) or Matrix (which is a local one).
 * Both Link and Matrix offer the same get and set functions.
 * Thus, the class "SU3" is a wrapper of different representation of SU3-matrices.
 */

// TODO using the typedefed "complex" type is not good style
template<class Type> class SU3
{
public:
	CUDA_HOST_DEVICE SU3( Type mat );
	CUDA_HOST_DEVICE SU3();
	Type mat;
	CUDA_HOST_DEVICE inline complex get( lat_group_dim_t i, lat_group_dim_t j );
	CUDA_HOST_DEVICE inline complex get(lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j);
	CUDA_HOST_DEVICE inline Quaternion<Real> getSubgroupQuaternion( lat_group_dim_t iSub, lat_group_dim_t jSub ); // TODO binding this class to class Quaternion is not good style -> make this a static function elsewhere
	CUDA_HOST_DEVICE inline void set( lat_group_dim_t i, lat_group_dim_t j, complex c);
	CUDA_HOST_DEVICE inline void set(lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j, complex c);
	CUDA_HOST_DEVICE inline SU3<Type>& operator+=( SU3<Type> ); // TODO overload for types like SU3<Link>::operator+=( SU3<Matrix> )
	template<bool reconstructThirdLine,class Type2> CUDA_HOST_DEVICE inline SU3<Type>& operator=( SU3<Type2> ); // TODO overload for types like SU3<Link>::operator+=( SU3<Matrix> )
	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type>& operator=( SU3<Type2> ); // TODO overload for types like SU3<Link>::operator+=( SU3<Matrix> )
	CUDA_HOST_DEVICE inline complex det();
	CUDA_HOST_DEVICE inline complex trace();
	CUDA_HOST_DEVICE inline void identity();
	CUDA_HOST_DEVICE inline void zero();
	CUDA_HOST_DEVICE inline void projectSU3();

	CUDA_HOST_DEVICE inline void reconstructThirdLine();

	CUDA_HOST_DEVICE inline void leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] );
	CUDA_HOST_DEVICE inline void rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] );

	CUDA_HOST_DEVICE inline void leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q );
	CUDA_HOST_DEVICE inline void rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Quaternion<Real> *q );

	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type> operator*( SU3<Type2> b  );
	template<class Type2> CUDA_HOST_DEVICE inline SU3<Type> operator+( SU3<Type2> b  );

	template<class SubgroupOperationClass> __device__ static inline void perSubgroup(SubgroupOperationClass t);
};

template<class Type> SU3<Type>::SU3( Type mat ) : mat(mat)
{
}

template<class Type> SU3<Type>::SU3()
{
}

template<class Type> complex SU3<Type>::get( lat_group_dim_t i, lat_group_dim_t j )
{
	return mat.get(i,j);
}

template<class Type> complex SU3<Type>::get( lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j )
{
	return mat.get( (i==0)?(iSub):(jSub), (j==1)?(iSub):(jSub) );	//Subgroup access TODO write about how this works
}

template<class Type> void SU3<Type>::set( lat_group_dim_t i, lat_group_dim_t j, complex c )
{
	return mat.set(i,j,c);
}

template<class Type> void SU3<Type>::set( lat_group_dim_t iSub, lat_group_dim_t jSub, lat_group_dim_t i, lat_group_dim_t j, complex c )
{
	return mat.set((i==0)?(iSub):(jSub), (j==1)?(iSub):(jSub) ,c);
}

/*
 * TODO SUBGROUP IS NOT PROJECTED TO SU2: this should be implied by its name
 */
template<class Type> Quaternion<Real> SU3<Type>::getSubgroupQuaternion( lat_group_dim_t iSub, lat_group_dim_t jSub )
{
	Quaternion<Real> q;
	complex temp;
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

	return q;
}


template<class Type> SU3<Type>& SU3<Type>::operator+=( SU3<Type> c )
{
	mat += c.mat;
	return *this;
}

template<class Type> template<class Type2> SU3<Type>& SU3<Type>::operator=( SU3<Type2> c )
{
	for( lat_group_dim_t i = 0; i < 2; i++ ) // TODO
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			mat.set(i,j,c.mat.get(i,j));
		}
	return *this;
}

template<class Type> template<bool reconstructThirdLine, class Type2> SU3<Type>& SU3<Type>::operator=( SU3<Type2> c )
{
	for( lat_group_dim_t i = 0; i < (reconstructThirdLine)?(2):(3); i++ )
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			mat.set(i,j,c.mat.get(i,j));
		}
	return *this;
}

template<class Type> template<class Type2> SU3<Type> SU3<Type>::operator*( SU3<Type2> b )
{
	SU3<Type> c;
	c.zero();
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for( lat_group_dim_t j = 0; j < 3; j++ )
		{
			complex temp(0,0);
			for( lat_group_dim_t k = 0; k < 3; k++ )
			{
				temp += get(i,k)*b.get(k,j);
			}
			c.set(i,j,temp);
		}
	}
	return c;
}

template<class Type> template<class Type2> SU3<Type> SU3<Type>::operator+( SU3<Type2> b )
{
	SU3<Type> c(*this);
	return c+=b;
}


template<class Type> complex SU3<Type>::det()
{
	complex c( 0, 0 );
	complex temp = get(0,0);
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

template<class Type> complex SU3<Type>::trace()
{
	return mat.trace();
}

template<class Type> void SU3<Type>::identity()
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for(lat_group_dim_t j = 0; j < 3; j++ )
		{
			if( i == j )
			{
				set( i,j, complex( 1, 0 ) );
			}
			else
			{
				set( i,j, complex( 0, 0 ) );
			}
		}
	}
}

template<class Type> void SU3<Type>::zero()
{
	for( lat_group_dim_t i = 0; i < 3; i++ )
	{
		for(lat_group_dim_t j = 0; j < 3; j++ )
		{
			set( i,j, complex(0, 0 ) );
		}
	}
}

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
	set( 2, 0, get(0,1).conj() * get(1,2).conj() - get(0,2).conj() * get(1,1).conj() );

	set( 2, 1, get(0,2).conj() * get(1,0).conj() - get(0,0).conj() * get(1,2).conj() );

	set( 2, 2, get(0,0).conj() * get(1,1).conj() - get(0,1).conj() * get(1,0).conj() );
}

template<class Type> void SU3<Type>::reconstructThirdLine()
{
	set( 2, 0, get(0,1).conj() * get(1,2).conj() - get(0,2).conj() * get(1,1).conj() );

	set( 2, 1, get(0,2).conj() * get(1,0).conj() - get(0,0).conj() * get(1,2).conj() );

	set( 2, 2, get(0,0).conj() * get(1,1).conj() - get(0,1).conj() * get(1,0).conj() );
}

/**
 * This actually works correctly.
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


template<class Type> void SU3<Type>::rightSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] )
{
	for( lat_group_dim_t k = 0; k < 3; k++ )
	{
//		ki = k*3+iSub;
		Complex<Real> KI = Complex<Real>(q[0],q[3]) * get(k,i);
		KI += Complex<Real>(-q[2],q[1]) * get(k,j);

		Complex<Real> KJ = Complex<Real>(q[2],q[1]) * get(k,i);
		KJ += Complex<Real>(q[0],-q[3]) * get(k,j);

		set( k, i , KI);
		set( k, j , KJ);
	}
}

template<class Type> void SU3<Type>::leftSubgroupMult( lat_group_dim_t i, lat_group_dim_t j, Real q[4] )
{
	for( lat_group_dim_t k = 0; k < 3; k++ )
	{
		Complex<Real> IK = Complex<Real>(q[0],q[3]) * get(i,k);
		IK += Complex<Real>(q[2],q[1])* get(j,k);

		Complex<Real> JK = Complex<Real>(-q[2], q[1]) * get(i,k);
		JK += Complex<Real>(q[0],-q[3]) * get(j,k);

		set( i, k , IK );
		set( j, k,  JK );
	}
}







template<class Type> template<class SubgroupOperationClass> void SU3<Type>::perSubgroup( SubgroupOperationClass t )
{
	t.subgroup(0,1);
	t.subgroup(0,2);
	t.subgroup(1,2);
}


#endif /* SU3_HXX_ */
