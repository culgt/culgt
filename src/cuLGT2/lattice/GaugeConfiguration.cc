#ifndef GAUGECONFIGURATION_CC
#define GAUGECONFIGURATION_CC

#include <exception>
#include <string>
#include "cuda.h"
#include "GaugeConfigurationHelper.h"

enum LocationType{ LocationHost, LocationDevice, LocationBoth };

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


template<typename PatternType> class GaugeConfiguration
{
public:
	// gathering lattice infos
	typedef typename PatternType::PARAMTYPE::TYPE T;
	static const int Ndim = PatternType::SITETYPE::Ndim;
	static const int LinkSize = PatternType::PARAMTYPE::SIZE;

	GaugeConfiguration( const int size[4] ) : UhostIsAllocated(false), UdeviceIsAllocated(false)
	{
		for( int i = 0; i < Ndim; i++ )
		{
			this->size[i] = size[i];
		}
		configurationSize = Ndim*getLatticeSize()*LinkSize;
	};


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


#endif


#include "gmock/gmock.h"
#include <iostream>

class GCSiteStub
{
public:
	static const int Ndim=4;
};

class GCParamStub
{
public:
	typedef double TYPE;
	static const int SIZE=18;
};

class GCPatternStub
{
public:
	typedef GCParamStub PARAMTYPE;
	typedef GCSiteStub SITETYPE;
};

class GaugeConfigurationDoubleFixedSize: public testing::Test {
public:
	int arraysize;
	int latticesize;
	static const int Nc = 3;
	static const int Ndim = 4;
	static const int Nt = 4;
	static const int Nx = 5;
	static const int Ny = 6;
	static const int Nz = 7;

	GaugeConfiguration<GCPatternStub>* gaugeconfig;

	void SetUp()
	{
		latticesize = Nt*Nx*Ny*Nz;
		arraysize = Nt*Nx*Ny*Nz*(Nc*Nc)*Ndim*2;
		int size[4] = {Nt,Nx,Ny,Nz};
		gaugeconfig = new GaugeConfiguration<GCPatternStub>( size );
	};

	void Teardown()
	{
		gaugeconfig->freeMemory();
		delete gaugeconfig;
		gaugeconfig = NULL;
	};
};


TEST_F( GaugeConfigurationDoubleFixedSize, LatticeSizesIsCorrectlySet )
{
	ASSERT_EQ( latticesize, gaugeconfig->getLatticeSize() );
}

TEST_F( GaugeConfigurationDoubleFixedSize, ConfigurationSizeIsCorrect )
{
	ASSERT_EQ( arraysize, gaugeconfig->getConfigurationSize() );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementFromHostThrowsExceptionIfNoMemoryIsAllocated )
{
	ASSERT_THROW( gaugeconfig->getElementFromHost( 3 ), MemoryException );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementFromHostDoesNotThrowExceptionIfMemoryIsAllocated )
{
	gaugeconfig->allocateMemory();

	ASSERT_NO_THROW( gaugeconfig->getElementFromHost( 3 ) );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementFromDeviceThrowsExceptionIfNoDeviceMemoryIsAllocated )
{
	gaugeconfig->allocateMemoryOnHost();
	ASSERT_THROW( gaugeconfig->getElementFromDevice( 3 ), MemoryException );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementFromDeviceDoesNotThrowExceptionIfMemoryIsAllocated )
{
	gaugeconfig->allocateMemoryOnDevice();
	ASSERT_NO_THROW( gaugeconfig->getElementFromDevice( 3 ) );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementNoThrowsAfterAllocateBoth )
{
	gaugeconfig->allocateMemory();
	ASSERT_NO_THROW( gaugeconfig->getElementFromHost( 3 ) );
	ASSERT_NO_THROW( gaugeconfig->getElementFromDevice( 3 ) );
}

TEST_F( GaugeConfigurationDoubleFixedSize, GetElementThrowsAfterFree )
{
	gaugeconfig->allocateMemory();
	gaugeconfig->freeMemory();
	ASSERT_THROW( gaugeconfig->getElementFromHost( 3 ), MemoryException );
	ASSERT_THROW( gaugeconfig->getElementFromDevice( 3 ), MemoryException );
}

TEST_F( GaugeConfigurationDoubleFixedSize, SetGetOnDeviceWorks )
{
	const int someIndex = 4;
	const double someValue = 1.5233;
	gaugeconfig->allocateMemoryOnDevice();
	gaugeconfig->setElementOnDevice( someIndex, someValue );

	ASSERT_DOUBLE_EQ( someValue, gaugeconfig->getElementFromDevice( someIndex ) );
}

TEST_F( GaugeConfigurationDoubleFixedSize, CopyToDeviceFailsIfNoMemoryIsAllocatedInHostAndDevice )
{
	gaugeconfig->allocateMemoryOnDevice();
	ASSERT_THROW( gaugeconfig->copyToDevice(), MemoryException );
}

TEST_F( GaugeConfigurationDoubleFixedSize, CopyConfigToDevice )
{
	const int someIndex = 4;
	const double someValue = 1.5233;
	gaugeconfig->allocateMemory();

	gaugeconfig->setElementOnHost( someIndex, someValue );
	gaugeconfig->copyToDevice();

	ASSERT_DOUBLE_EQ( someValue, gaugeconfig->getElementFromDevice( someIndex ) );
}

TEST_F( GaugeConfigurationDoubleFixedSize, CopyConfigToHost )
{
	const int someIndex = 46;
	const double someValue = 2.5345;
	gaugeconfig->allocateMemory();

	gaugeconfig->setElementOnDevice( someIndex, someValue );
	gaugeconfig->copyToHost();

	ASSERT_DOUBLE_EQ( someValue, gaugeconfig->getElementFromHost( someIndex ) );
}

