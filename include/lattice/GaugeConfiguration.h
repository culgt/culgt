#ifndef GAUGECONFIGURATION_CC
#define GAUGECONFIGURATION_CC

#include <cuda.h>
#include <cuda_runtime.h>
#include <exception>
#include <string>
#include "GaugeConfigurationCudaHelper.h"
#include "LocalLink.h"
#include "GlobalLink.h"
#include "filetypes/LinkFile.h"
#include "LatticeDimension.h"

namespace culgt {

class MemoryException: public std::exception
{
public:
	~MemoryException() throw() {};
	MemoryException( std::string msg ) : msg(msg){};
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
private:
	std::string msg;
};

template<typename PatternType> class GaugeConfiguration
{
public:
	// gathering lattice infos
	typedef PatternType PATTERNTYPE;
	typedef typename PatternType::PARAMTYPE::TYPE T;
	typedef typename PatternType::PARAMTYPE PARAMTYPE;
	typedef typename PatternType::SITETYPE SITETYPE;
	static const int Ndim = PatternType::SITETYPE::NDIM;
	static const int LinkSize = PatternType::PARAMTYPE::SIZE;

	GaugeConfiguration( const LatticeDimension<Ndim> dim ) : Uhost(NULL), Udevice(NULL), UhostIsAllocated(false), UdeviceIsAllocated(false), dim(dim)
	{
		configurationSize = Ndim*dim.getSize()*LinkSize;
	};

	~GaugeConfiguration()
	{
		freeMemory();
	}


	void allocateMemory()
	{
		allocateMemoryOnHost();
#ifdef __CUDACC__
		allocateMemoryOnDevice();
#endif
	}

	void allocateMemoryOnHost()
	{
		if( !UhostIsAllocated )
		{
			Uhost = new T[configurationSize];
			UhostIsAllocated = true;
		}
	}

	void freeMemory()
	{
		freeMemoryOnHost();
#ifdef __CUDACC__
		freeMemoryOnDevice();
#endif
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


	LatticeDimension<Ndim> getLatticeDimension() const
	{
		return dim;
	}

	int getConfigurationSize() const
	{
		return configurationSize;
	}

	void loadFile( LinkFile<PatternType>& linkfile )
	{
		linkfile.load( Uhost );
	}


	T* getHostPointer() const
	{
		return Uhost;
	}

#ifdef __CUDACC__
	void allocateMemoryOnDevice()
	{
		if( !UdeviceIsAllocated )
		{
			GaugeConfigurationCudaHelper<T>::allocateMemory( &Udevice, configurationSize );
			UdeviceIsAllocated = true;
		}
	}

	void copyToHost()
	{
		checkMemoryIsAllocated();
		GaugeConfigurationCudaHelper<T>::copyToHost( Uhost, Udevice, configurationSize );
	}

	void freeMemoryOnDevice()
	{
		if( UdeviceIsAllocated )
		{
			GaugeConfigurationCudaHelper<T>::freeMemory( Udevice );
			Udevice = NULL;
			UdeviceIsAllocated = false;
		}
	}

	void copyToDevice()
	{
		checkMemoryIsAllocated();
		GaugeConfigurationCudaHelper<T>::copyToDevice( Udevice, Uhost, configurationSize );
	}

	T* getDevicePointer() const
	{
		return Udevice;
	}

	T* getDevicePointer( lat_coord_t timeslice )
	{
		BOOST_MPL_ASSERT_RELATION( SITETYPE::PARITYTYPE, ==, TIMESLICE_SPLIT );
		lat_array_index_t configurationSizeTimeslice = Ndim*dim.getSizeTimeslice()*LinkSize;
		return &Udevice[timeslice*configurationSizeTimeslice];
	}

	void setDevicePointer( T* U, bool assumeAllocated = false )
	{
		Udevice = U;
		UdeviceIsAllocated = assumeAllocated;
	}

	void takeCopyDeviceToDevice( GaugeConfiguration<PatternType>& src )
	{
		GaugeConfigurationCudaHelper<T>::copyDeviceToDevice( Udevice, src.getDevicePointer(), configurationSize );
	}

	T getElementFromDevice( int i ) const
	{
		if( UdeviceIsAllocated )
			return GaugeConfigurationCudaHelper<T>::getElement( Udevice, i );
		else
			throw MemoryException( "No device memory allocated!" );
	}

	void setElementOnDevice( int i, const T val )
	{
		if( UdeviceIsAllocated )
			GaugeConfigurationCudaHelper<T>::setElement( Udevice, i, val );
		else
			throw MemoryException( "No device memory allocated!" );
	}

	void setLinkOnDevice( SITETYPE site, int mu, const LocalLink<PARAMTYPE> link )
	{
		if( UdeviceIsAllocated )
		{
			GaugeConfigurationCudaHelper<T>::template setLink<PatternType>( Udevice, site, mu, link );
		}
		else
			throw MemoryException( "No device memory allocated!" );
	}

	LocalLink<PARAMTYPE> getLinkFromDevice( SITETYPE site, int mu ) const
	{
		if( UdeviceIsAllocated )
		{
			return GaugeConfigurationCudaHelper<T>::template getLink<PatternType,LocalLink<PARAMTYPE> >( Udevice, site, mu );
		}
		else
		{
			throw MemoryException( "No device memory allocated!" );
		}
	}

	template<typename RNG> void setHotOnDevice( int rngSeed, int rngCounter )
	{
		GaugeConfigurationCudaHelper<T>::template setHot<PatternType,RNG>(Udevice, dim, rngSeed, rngCounter );
	}
	void setColdOnDevice()
	{
		GaugeConfigurationCudaHelper<T>::template setCold<PatternType>(Udevice, dim );
	}
#endif

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
	LatticeDimension<Ndim> dim;
	lat_array_index_t configurationSize;
};


}

#endif



