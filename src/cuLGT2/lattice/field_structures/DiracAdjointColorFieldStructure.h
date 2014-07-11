/**
 *
 * Memory layout is (slowest running index first): dirac index(mu), color index (a), site index
 *
 *  Created on: Jun 23, 2014
 *      Author: vogt
 */


#ifndef DIRACADJOINTCOLORFIELDSTRUCTURE_H_
#define DIRACADJOINTCOLORFIELDSTRUCTURE_H_

#include "../../common/culgt_typedefs.h"
#include "FieldStructure.h"
#include "../../cudacommon/cuda_host_device.h"

namespace culgt
{

template<typename SiteType, int Nc> class DiracAdjointColorFieldStructure: public FieldStructure
{
public:
	typedef SiteType SITETYPE;

	CUDA_HOST_DEVICE DiracAdjointColorFieldStructure( int mu, int a, SiteType site ) : mu(mu), a(a), site( site )
	{
	}

	CUDA_HOST_DEVICE DiracAdjointColorFieldStructure( SiteType site ) :  mu(0), a(0), site( site )
	{
	}

	CUDA_HOST_DEVICE lat_array_index_t getIndex() const
	{
		return (mu*(Nc*Nc-1)+a)*site.getSize()+site.getIndex();
	}

	CUDA_HOST_DEVICE lat_array_index_t getSize() const
	{
		return (Nc*Nc-1)*SiteType::Ndim*site.getSize();
	}

	CUDA_HOST_DEVICE SiteType& getSite()
	{
		return site;
	}

	CUDA_HOST_DEVICE int getA() const
	{
		return a;
	}

	CUDA_HOST_DEVICE void setA(int a)
	{
		this->a = a;
	}

	CUDA_HOST_DEVICE int getMu() const
	{
		return mu;
	}

	CUDA_HOST_DEVICE void setMu(int mu)
	{
		this->mu = mu;
	}

private:
	int mu;
	int a;
	SiteType site;
};

}

#endif
