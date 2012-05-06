/*
 * GpuTimeslicePattern.hxx
 *
 *  Created on: Apr 13, 2012
 *      Author: vogt
 */

#ifndef GPULANDAUPATTERN_HXX_
#define GPULANDAUPATTERN_HXX_

#include "../../util/cuda/cuda_host_device.h"
#include "../SiteCoord.hxx"
#include <assert.h>

template<class Site, int T_Ndim, int T_Nc> class GpuLandauPattern
{
public:
	static const int Nc = T_Nc;
	static const int Ndim = T_Ndim;
	CUDA_HOST_DEVICE static inline int getSiteIndex( Site s );
	CUDA_HOST_DEVICE static inline int getLinkIndex( Site s, int mu );
	CUDA_HOST_DEVICE static inline int getIndex( int linkIndex, int i, int j, int c );
	CUDA_HOST_DEVICE static inline int getIndex( Site s, int mu, int i, int j, int c );
};

template<class Site, int T_Ndim, int T_Nc> int GpuLandauPattern<Site, T_Ndim,T_Nc>::getSiteIndex( Site s )
{
	return  s.getLatticeIndex();
}

template<class Site, int T_Ndim, int T_Nc> int GpuLandauPattern<Site, T_Ndim,T_Nc>::getLinkIndex( Site s, int mu )
{
	int muSize = T_Nc*T_Nc*2*s.getLatticeSize();
	return s.getLatticeIndex()+mu*muSize;
}

template<class Site, int T_Ndim, int T_Nc> int GpuLandauPattern<Site, T_Ndim,T_Nc>::getIndex( Site s, int mu, int i, int j, int c )
{
	return s.getLatticeIndex() + s.getLatticeSize()*( c + 2 * ( j + T_Nc *( i + T_Nc * mu ) ) );
}

#endif /* GPUTIMESLICEPATTERN_HXX_ */
