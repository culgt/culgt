/**
 * Defines integer types for different purposes, like lattice coordinates, lattice index, array index (for the global link array)
 * This is included in every file!
 *
 * If you want to change the settings for a specific program include a copy of this file with your settings before everything else is included!
 */

#ifndef LATTICE_TYPEDEFS_H_
#define LATTICE_TYPEDEFS_H_

#include "../cudacommon/cuda_host_device.h"
#include "cuda_functions.h"

#if !defined(__VECTOR_TYPES_H__)
	struct float4 {
		float x;
		float y;
		float z;
		float w;
	};
	struct double4 {
		double x;
		double y;
		double z;
		double w;
	};
	struct float2 {
		float x;
		float y;
	};
	struct double2 {
		double x;
		double y;
	};
#endif

#if __CUDA_ARCH__ >= 350
#if (defined(_MSC_VER) && defined(_WIN64)) || defined(__LP64__)
#define __LDG_PTR   "l"
#else
#define __LDG_PTR   "r"
#endif

static __device__ __inline__ double4 __ldg(const double4 *ptr) { const double2* ptr2 = reinterpret_cast<const double2*>(ptr); double4 ret; asm volatile ("ld.global.nc.v2.f64 {%0,%1}, [%2];"  : "=d"(ret.x), "=d"(ret.y) : __LDG_PTR (ptr2)); ptr2++; asm volatile ("ld.global.nc.v2.f64 {%0,%1}, [%2];"  : "=d"(ret.z), "=d"(ret.w) : __LDG_PTR (ptr2)); return ret; }

#undef __LDG_PTR
#endif

template<typename T> struct Real4
{
};

template<> struct Real4<float>
{
	typedef float4 VECTORTYPE;
};

template<> struct Real4<double>
{
	typedef double4 VECTORTYPE;
};

template<typename T> struct Real2
{
};

template<> struct Real2<float>
{
	typedef float2 VECTORTYPE;
};

template<> struct Real2<double>
{
	typedef double2 VECTORTYPE;
};




inline bool operator==(float4 const& lhs, float4 const& rhs)
{
	return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}
inline bool operator!=(float4 const& lhs, float4 const& rhs)
{
	return !(lhs == rhs);
}
inline bool operator==(double4 const& lhs, double4 const& rhs)
{
	return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}
inline bool operator!=(double4 const& lhs, double4 const& rhs)
{
	return !(lhs == rhs);
}


// type for lattic dimension, i.e. 3 or 4, ...
//  and for group dimensions
typedef int lat_dim_t;

// type for things like parity (maybe there could be a benefit in a compiler when replacing this with bool, but probably not)
typedef int lat_bool_t;

// type for group indices, i.e. the index in the parameterization (for float18 this would be 0..17)
typedef int lat_group_index_t;

// datatype for lattice coordinates, i.e. the x0, x1, ... coordinates, as wells as the corresponding lattice sizes
//typedef short lat_coord_t;
typedef int lat_coord_t;

// type for the lattice index, i.e. the index built from the coordinates
typedef int lat_index_t;

// type for index of the global link array (for really large arrays this might need long int)
typedef int lat_array_index_t;
//typedef unsigned int lat_array_index_t;

#endif /* LATTICE_TYPEDEFS_H_ */
