#ifndef GAUGECONFIGURATION_CC
#define GAUGECONFIGURATION_CC

#include <cuda.h>
#include <cuda_runtime.h>
#include <exception>
#include <string>
#include "GaugeConfigurationHelper.h"
#include "LocalLink.h"
#include "GlobalLink.h"
#include "LinkFile.h"

namespace culgt {

class MemoryException: public std::exception
{
public:
	MemoryException( std::string msg ) : msg(msg){};
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
private:
	std::string msg;
};

/**
 * TODO Maybe it would be good to have a method where you can pass in GaugeConfiguration with other pattern.
 */

template<typename PatternType> class GaugeConfiguration
{
public:
	// gathering lattice infos
	typedef typename PatternType::PARAMTYPE::TYPE T;
	typedef typename PatternType::PARAMTYPE PARAMTYPE;
	typedef typename PatternType::SITETYPE SITETYPE;
	static const int Ndim = PatternType::SITETYPE::Ndim;
	static const int LinkSize = PatternType::PARAMTYPE::SIZE;

	GaugeConfiguration( const int size[Ndim] ) : UhostIsAllocated(false), UdeviceIsAllocated(false)
	{
		for( int i = 0; i < Ndim; i++ )
		{
			this->size[i] = size[i];
		}
		configurationSize = Ndim*getLatticeSize()*LinkSize;
	};

	~GaugeConfiguration()
	{
		freeMemory();
	}


	void allocateMemory()
	{
		allocateMemoryOnHost();
		allocateMemoryOnDevice();
	}

	void allocateMemoryOnHost()
	{
		if( !UhostIsAllocated )
		{
			Uhost = new T[configurationSize];
			UhostIsAllocated = true;
		}
	}

	void allocateMemoryOnDevice()
	{
		if( !UdeviceIsAllocated )
		{
			GaugeConfigurationHelper<T>::allocateMemory( &Udevice, configurationSize );
			UdeviceIsAllocated = true;
		}
	}

	void freeMemory()
	{
		freeMemoryOnHost();
		freeMemoryOnDevice();
	}

	void freeMemoryOnHost()
	{
		if( UhostIsAllocated )
		{
			delete[] Uhost;
			Uhost = NULL;
			UhostIsAllocated = false;
		}
	}

	void freeMemoryOnDevice()
	{
		if( UdeviceIsAllocated )
		{
			GaugeConfigurationHelper<T>::freeMemory( Udevice );
			Udevice = NULL;
			UdeviceIsAllocated = false;
		}
	}

	void copyToDevice()
	{
		checkMemoryIsAllocated();
		GaugeConfigurationHelper<T>::copyToDevice( Udevice, Uhost, configurationSize );
	}

	void copyToHost()
	{
		checkMemoryIsAllocated();
		GaugeConfigurationHelper<T>::copyToHost( Uhost, Udevice, configurationSize );
	}

	T getElementFromHost( int i ) const
	{
		if( UhostIsAllocated )
			return Uhost[i];
		else
			throw MemoryException( "No host memory allocated!" );
	}

	void setElementOnHost( int i, const T val )
	{
		if( UhostIsAllocated )
			Uhost[i] = val;
		else
			throw MemoryException( "No host memory allocated!" );
	}

	void setLinkOnHost( SITETYPE site, int mu, const LocalLink<PARAMTYPE> link )
	{
		if( UhostIsAllocated )
		{
			GlobalLink<PatternType> globalLink( Uhost, site, mu );
			globalLink = link;
		}
		else
			throw MemoryException( "No host memory allocated!" );
	}

	LocalLink<PARAMTYPE> getLinkFromHost( SITETYPE site, int mu ) const
	{
		if( UhostIsAllocated )
		{
			LocalLink<PARAMTYPE> link;
			GlobalLink<PatternType> globalLink( Uhost, site, mu );
			link = globalLink;
			return link;
		}
		else
			throw MemoryException( "No host memory allocated!" );
	}

	T getElementFromDevice( int i ) const
	{
		if( UdeviceIsAllocated )
			return GaugeConfigurationHelper<T>::getElement( Udevice, i );
		else
			throw MemoryException( "No device memory allocated!" );
	}

	void setElementOnDevice( int i, const T val )
	{
		if( UdeviceIsAllocated )
			GaugeConfigurationHelper<T>::setElement( Udevice, i, val );
		else
			throw MemoryException( "No device memory allocated!" );
	}

	void setLinkOnDevice( SITETYPE site, int mu, const LocalLink<PARAMTYPE> link )
	{
		if( UdeviceIsAllocated )
		{
			// TODO this is an ugly construct...
			GaugeConfigurationHelper<T>::template setLink<PatternType>( Udevice, site, mu, link );
		}
		else
			throw MemoryException( "No device memory allocated!" );
	}

	LocalLink<PARAMTYPE> getLinkFromDevice( SITETYPE site, int mu ) const
	{
		if( UdeviceIsAllocated )
		{
			return GaugeConfigurationHelper<T>::template getLink<PatternType,LocalLink<PARAMTYPE> >( Udevice, site, mu );;
		}
		else
			throw MemoryException( "No device memory allocated!" );
	}

	int getLatticeSize() const
	{
		int latsize = 1;
		for( int i = 0; i < Ndim; i++ )
		{
			latsize *= size[i];
		}
		return latsize;
	}

	int getConfigurationSize() const
	{
		return configurationSize;
	}

	void loadFile( LinkFile<PatternType>& linkfile )
	{
		linkfile.load( Uhost );
	}

private:
	void checkMemoryIsAllocated()
	{
		if( !UhostIsAllocated )
			throw MemoryException( "No host memory allocated!" );
		else if( !UdeviceIsAllocated )
			throw MemoryException( "No device memory allocated!" );
	}

	T* Uhost;
	T* Udevice;
	bool UhostIsAllocated;
	bool UdeviceIsAllocated;
	int size[Ndim];
	int configurationSize;
};


}

#endif



