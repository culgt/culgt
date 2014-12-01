/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef GLOBALLINK_H_
#define GLOBALLINK_H_

#include "../cudacommon/cuda_host_device.h"
#include "../common/culgt_typedefs.h"
#include "ParameterizationMediator.h"
#include "Link.h"
#include "LatticeDimension.h"
#include <boost/static_assert.hpp>

#ifdef __CUDACC__
#include "TextureManager.h"

#define COPY_GLOBALLINKTYPE( SRC, DEST, TEXTUREID )\
	typedef GlobalLink<typename SRC::PATTERNTYPE, SRC::USETEXTURE, TEXTUREID> DEST

#endif

namespace culgt
{
template<typename ConfigurationPattern, bool UseTexture=false, int TextureID = 0> class GlobalLink//: Link<typename ConfigurationPattern::PARAMTYPE>
{
public:

	/**
	 * Layout of the Link (for example 18 real parameters, 9 complex parameters, or 12 real, ...)
	 */
	typedef typename ConfigurationPattern::PARAMTYPE PARAMTYPE;
	typedef ConfigurationPattern PATTERNTYPE;
	typedef ConfigurationPattern CONFIGURATIONPATTERN;
	static const bool USETEXTURE=UseTexture;


#ifdef __CUDACC__
	BOOST_STATIC_ASSERT_MSG( TextureID < TextureManager<typename PARAMTYPE::TYPE>::MAX_TEXTURES, "Max. texture number exceeded, please add more textures in TextureManager.h");
	static cudaError_t bindTexture( typename PARAMTYPE::TYPE* pointerToStore, lat_array_index_t arraySize )
	{
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<typename TextureManager<typename PARAMTYPE::TYPE>::Type>();
		size_t offset;
		return cudaBindTexture(&offset, TextureManager<typename PARAMTYPE::TYPE>::getTextureReference( TextureID ), pointerToStore, channelDesc, arraySize*sizeof(typename PARAMTYPE::TYPE) );
	}
	static void unbindTexture()
	{
		cudaUnbindTexture( TextureManager<typename PARAMTYPE::TYPE>::getTextureReference( TextureID ) );
	}
#endif

	CUDA_HOST_DEVICE inline ~GlobalLink()
	{
	};

	CUDA_HOST_DEVICE inline GlobalLink( typename PARAMTYPE::TYPE* pointerToStore, const typename ConfigurationPattern::SITETYPE site, lat_dim_t mu ) : pointerToStore(pointerToStore), site(site), mu(mu) {} ;

	/**
	 * Returns the i-th entry of the ParamType specific array
	 *
	 *  I want something like: store[accesspattern(site,mu,i)]
	 *	where
	 *	- site is an object referring to a lattice site
	 *	- mu is the direction of the link
	 *	- i is the index that belongs to a parameterization
	 *	  (for example with parameterization "float12": i=0 is the real part of the 1/1, i = 1 the imaginary part, i=2 real of 1/2 and so on)
	 * @param i
	 * @return
	 */
	CUDA_HOST_DEVICE inline typename PARAMTYPE::TYPE get( lat_group_index_t i ) const
	{
#ifdef __CUDA_ARCH__
		if( UseTexture )
		{
			return Tex1DFetcher<typename PARAMTYPE::TYPE>::fetch( TextureManager<typename PARAMTYPE::TYPE>::getTexture( TextureID ), ConfigurationPattern::getIndex(site,mu,i) );;
		}
		else
#endif
		{
			return pointerToStore[ConfigurationPattern::getIndex(site,mu,i)];

		}
	}
	/**
	 * Sets the i-th entry of the ParamType specific array
	 * @param i
	 * @param val
	 */
	CUDA_HOST_DEVICE inline void set( lat_group_index_t i, typename PARAMTYPE::TYPE val)
	{
		pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = val;
	}

	CUDA_HOST_DEVICE inline void zero()
	{
		for( lat_group_index_t i = 0; i < PARAMTYPE::SIZE; i++ )
		{
			pointerToStore[ConfigurationPattern::getIndex(site,mu,i)] = (typename PARAMTYPE::TYPE) 0.0;
		}
	}

	/**
	 * Assignment operator for same type just copies all elements.
	 * @param arg
	 * @return
	 */
	CUDA_HOST_DEVICE inline GlobalLink<ConfigurationPattern,USETEXTURE>& operator=( const GlobalLink<ConfigurationPattern,USETEXTURE>& arg )
	{
		for( lat_group_index_t i = 0; i < ConfigurationPattern::PARAMTYPE::SIZE; i++ )
		{
			set( i, arg.get(i) );
		}
		return *this;
	}

//	CUDA_HOST_DEVICE inline GlobalLink<ConfigurationPattern,USETEXTURE>& operator=( const LocalLink<typename ConfigurationPattern::PARAMTYPE>& arg )
//	{
//		for( lat_group_index_t i = 0; i < ConfigurationPattern::PARAMTYPE::SIZE; i++ )
//		{
//			set( i, arg.get(i) );
//		}
//		return *this;
//	}

	static inline lat_array_index_t getArraySize( LatticeDimension<ConfigurationPattern::SITETYPE::Ndim> dim )
	{
		return dim.getSize()*ConfigurationPattern::SITETYPE::Ndim*PARAMTYPE::SIZE;
	}

	/**
	 * Assignment operator needs to call Mediator for different ParamTypes and/or LinkTypes
	 * @param arg
	 * @return
	 */
	template<typename LinkType> inline CUDA_HOST_DEVICE GlobalLink<ConfigurationPattern,USETEXTURE,TextureID>& operator=( const LinkType arg )
	{
		ParameterizationMediator<typename ConfigurationPattern::PARAMTYPE,typename LinkType::PARAMTYPE,GlobalLink<ConfigurationPattern,USETEXTURE,TextureID>, LinkType >::assign( *this, arg );
		return *this;
	}

	typename PARAMTYPE::TYPE* pointerToStore;
	const typename ConfigurationPattern::SITETYPE site;
	lat_dim_t mu;
};


} /* namespace culgt */

#endif /* GLOBALLINK_H_ */



