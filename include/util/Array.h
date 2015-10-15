/*
 */

#ifndef ARRAY_H_
#define ARRAY_H_

namespace culgt
{

template<typename T, int Size> class Array
{
private:
	T data[Size];
public:
	CUDA_HOST_DEVICE inline Array()
	{
	}
	CUDA_HOST_DEVICE inline Array( T val )
	{
		for( int i = 0; i < Size; ++i ) data[i] = val;
	}
	CUDA_HOST_DEVICE inline T operator[]( int i ) const
	{
		return data[i];
	}
	CUDA_HOST_DEVICE inline T& operator[]( int i )
	{
		return data[i];
	}

	CUDA_HOST_DEVICE inline Array& operator=( const Array& src )
	{
		for( int i = 0; i < Size; ++i )
		{
			data[i] = src[i];
		}
		return *this;
	}

};

}



#endif /* ARRAY_H_ */
