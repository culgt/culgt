#include "gmock/gmock.h"
#include "lattice/GaugeConfiguration.h"
#include "lattice/filetypes/LinkFile.h"
#include "testhelper_pattern_stub.h"
using namespace culgt;



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
	static const int size[4];

	GaugeConfiguration<PatternStub<double> >* gaugeconfig;

	void SetUp()
	{
		latticesize = Nt*Nx*Ny*Nz;
		arraysize = Nt*Nx*Ny*Nz*(Nc*Nc)*Ndim*2;

		gaugeconfig = new GaugeConfiguration<PatternStub<double> >( size );
	};

	void Teardown()
	{
		gaugeconfig->freeMemory();
		delete gaugeconfig;
		gaugeconfig = NULL;
	};
};

const int GaugeConfigurationDoubleFixedSize::size[4] = {Nt,Nx,Ny,Nz};

TEST_F( GaugeConfigurationDoubleFixedSize, LatticeSizesIsCorrectlySet )
{
	ASSERT_EQ( latticesize, gaugeconfig->getLatticeDimension().getSize() );
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

TEST_F( GaugeConfigurationDoubleFixedSize, CopyToDeviceFailsIfNoMemoryIsAllocatedInHost )
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



class LinkFileMock: public LinkFile<PatternStub<float> >
{
public:
	MOCK_METHOD0(loadImplementation, void() );
};

class GaugeConfigurationFileLoad: public testing::Test
{
public:
	static const int size[4];
	GaugeConfiguration<PatternStub<float> > gaugeconfig;
	LinkFileMock linkfile;

	GaugeConfigurationFileLoad() : gaugeconfig(size)
	{
	}

	void SetUp()
	{
		linkfile.setFilename( "/dev/random" );
	}
};

const int GaugeConfigurationFileLoad::size[4] = {4,4,4,4};

TEST_F( GaugeConfigurationFileLoad, CallsLoadImplementation )
{
	EXPECT_CALL( linkfile, loadImplementation() );

	gaugeconfig.loadFile( linkfile );
}
