/**
 * TextureManager.h
 *
 * Texture references need to be global variables. We define them here and use template overloads to pick the correct reference.
 * Tex1DFetcher overloads the tex1Dfetch to allow double and double4 texture loads
 *
 * Allows up to 2 texture references of the same type. If we need more at some point: do some Meta programming magic...
 */

#ifndef TEXTUREMANAGER_H_
#define TEXTUREMANAGER_H_

#include "../cudacommon/cuda_host_device.h"

namespace culgt
{

texture<float4> float4texture_0;
texture<float4> float4texture_1;
texture<int4> double4texture_0;
texture<int4> double4texture_1;
texture<float> floattexture_0;
texture<float> floattexture_1;
texture<int2> doubletexture_0;
texture<int2> doubletexture_1;

template<typename T> class TextureManager
{
};

template<> class TextureManager<float4>
{
public:
	static const int MAX_TEXTURES = 2;
	typedef texture<float4> TextureType;
	typedef float4 Type;
	static inline CUDA_HOST_DEVICE TextureType getTexture(int textureID = 0)
	{
		if( textureID == 0)
		{
			return float4texture_0;
		}
		else
		{
			return float4texture_1;
		}
	}
	static inline __host__ TextureType& getTextureReference(int textureID = 0)
	{
		if( textureID == 0)
		{
			return float4texture_0;
		}
		else
		{
			return float4texture_1;
		}
	}
};

template<> class TextureManager<double4>
{
public:
	static const int MAX_TEXTURES = 2;
	typedef texture<int4> TextureType;
	typedef int4 Type;
	static inline CUDA_HOST_DEVICE TextureType getTexture(int textureID = 0)
	{
		if( textureID == 0)
		{
			return double4texture_0;
		}
		else
		{
			return double4texture_1;
		}
	}
	static inline __host__ TextureType& getTextureReference(int textureID = 0)
	{
		if( textureID == 0)
		{
			return double4texture_0;
		}
		else
		{
			return double4texture_1;
		}
	}
};

template<> class TextureManager<float>
{
public:
	static const int MAX_TEXTURES = 2;
	typedef texture<float> TextureType;
	typedef float Type;
	static inline CUDA_HOST_DEVICE TextureType getTexture(int textureID = 0)
	{
		if( textureID == 0)
		{
			return floattexture_0;
		}
		else
		{
			return floattexture_1;
		}
	}
	static inline __host__ TextureType& getTextureReference(int textureID = 0)
	{
		if( textureID == 0)
		{
			return floattexture_0;
		}
		else
		{
			return floattexture_1;
		}
	}
};

template<> class TextureManager<double>
{
public:
	static const int MAX_TEXTURES = 2;
	typedef texture<int2> TextureType;
	typedef int2 Type;
	static inline CUDA_HOST_DEVICE TextureType getTexture(int textureID = 0)
	{
		if( textureID == 0)
		{
			return doubletexture_0;
		}
		else
		{
			return doubletexture_1;
		}
	}
	static inline __host__ TextureType& getTextureReference(int textureID = 0)
	{
		if( textureID == 0)
		{
			return doubletexture_0;
		}
		else
		{
			return doubletexture_1;
		}
	}
};


template<typename T> class Tex1DFetcher
{
public:
	static inline CUDA_HOST_DEVICE T fetch( texture<T> tex, lat_array_index_t index )
	{
		return tex1Dfetch( tex, index );
	}
};

template<> class Tex1DFetcher<double4>
{
public:
	static inline __device__ double4 fetch( texture<int4> tex, lat_array_index_t index )
	{
		double4 result;
		int4 temp;
		temp = tex1Dfetch( tex, index*2 );
		result.x = __hiloint2double(temp.y, temp.x);
		result.y = __hiloint2double(temp.w, temp.z);
		temp = tex1Dfetch( tex, index*2+1 );
		result.z = __hiloint2double(temp.y, temp.x);
		result.w = __hiloint2double(temp.w, temp.z);
		return result;
	}
};

template<> class Tex1DFetcher<double>
{
public:
	static inline __device__ double fetch( texture<int2> tex, lat_array_index_t index )
	{
		double result;
		int2 temp;
		temp = tex1Dfetch( tex, index );
		result = __hiloint2double(temp.y, temp.x);
		return result;
	}
};

}

#endif /* TEXTUREMANAGER_H_ */
