/**
 * Quaternion representation for a SU2 (or proportional) element.
 *
 * TODO:
 *  - implement all reasonable operations: *,+,...
 *  - This may be a storage class of the not yet implemented "SU2" frontend class. Wenn tackling the question of where to
 *    place the mathematical operations (in frontend or storage classes) consider the comment of the trace() function.
 *
 * For an introduction to the relation of SU2 and quaternion see for example Gattringer & Lang (2010), p. 81.
 * The first row of a SU2-matrix is given by ( x[0]+i*x[3]    x[2]+i*x[1] ).
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-04-23
 */

#ifndef QUATERNION_HXX_
#define QUATERNION_HXX_

template<class T> class Quaternion
{
public:
	CUDA_HOST_DEVICE inline Quaternion();
	CUDA_HOST_DEVICE inline virtual ~Quaternion();
	CUDA_HOST_DEVICE inline T& operator[]( int i );
	CUDA_HOST_DEVICE inline Complex<T> get( int i, int j );
	CUDA_HOST_DEVICE inline void set( int i, int j, Complex<T> c );
	CUDA_HOST_DEVICE inline Quaternion<T>& hermitian();
	CUDA_HOST_DEVICE inline Complex<T> det(); //is the implementation correct???
	CUDA_HOST_DEVICE inline Quaternion<T>& project_SU2();
private:
	T quat[4];
};

template<class T> Quaternion<T>::Quaternion()
{
}

template<class T> Quaternion<T>::~Quaternion()
{
}

/**
 * Return element of quaternion representation.
 * @parameter index i=0,1,2,3
 * @return element
 */
template<class T> T& Quaternion<T>::operator[]( int i )
{
	return quat[i];
}


/**
 * Returns the element (i,j) if the quaternion is considered as a storage type for SU2.
 * For example get(0,0) returns the complex value x[0]+i*x[3], get(1,1) returns x[0]-i*x[3].
 *
 * TODO element access by this pattern through class "SU2" is not a good idea. -> to many calls to get() for a trace() for example...
 * specialization of SU2 for quaternion type is needed here!!!
 */
template<class T> Complex<T> Quaternion<T>::get( int i, int j )
{
	switch( i )
	{
	case 0:
		switch( j )
		{
		case 0:
			return Complex<T>( quat[0], quat[3] );
		case 1:
			return Complex<T>( quat[2], quat[1] );
		}
	case 1:
		switch( j )
		{
		case 0:
			return Complex<T>( -quat[2], quat[1] );
		case 1:
			return Complex<T>( quat[0], -quat[3] );
		}
	}
	return Complex<T>(0,0);
}

/**
 * TODO see get(); here it is is even more obvious: if we set a SU2 matrix to some value,
 * we do the set in Quaternion twice if we are not aware of the underlying storage class
 */
template<class T> void Quaternion<T>::set( int i, int j, Complex<T> c )
{
	switch( i )
	{
	case 0:
		switch( j )
		{
		case 0:
			quat[0] = c.x;
			quat[3] = c.y;
			break;
		case 1:
			quat[2] = c.x;
			quat[1] = c.y;
			break;
		}
	case 1:
		switch( j )
		{
		case 0:
			quat[2] = -c.x;
			quat[1] = c.y;
			break;
		case 1:
			quat[0] = c.x;
			quat[3] = -c.y;
		}
	}
}

/**
 * Computes the hermitian conjugate.
 */
template<class T> Quaternion<T>& Quaternion<T>::hermitian()
{
	quat[1] = -quat[1];
	quat[2] = -quat[2];
	quat[3] = -quat[3];
	return *this;
}

/**
 * Computes the det.
 * TODO implement the imaginary part!
 */
template<class T> Complex<T> Quaternion<T>::det()
{
	Complex<T> c( quat[0]*quat[0]-quat[3]*quat[3]+quat[2]*quat[2]-quat[1]*quat[1], 0 ); // TODO im
	return c;
}

/**
 * normalizes to det=1
 */
template<class T> Quaternion<T>& Quaternion<T>::project_SU2()
{
	T c = rsqrt(quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]);
	quat[0] *= c;
	quat[1] *= c;
	quat[2] *= c;
	quat[3] *= c;
	return *this;
}

#endif /* QUATERNION_HXX_ */
