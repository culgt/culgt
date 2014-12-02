/**
 * test_DeviceReader.cu
 *
 *  Created on: Mar 14, 2014
 *      Author: vogt
 */


#include "gmock/gmock.h"
#include "cudacommon/DeviceCommunicator.h"
#include "cudacommon/cuda_error.h"

using namespace culgt;
using namespace ::testing;

class ADeviceCommunicator: public Test
{
public:
	float* cudaMemory;
//	static constexpr float someValue = 1.523;
//	static constexpr int someIndex = 4;
	static const float someValue = 1.523;
	static const int someIndex = 4;

	void SetUp()
	{
		CUDA_SAFE_CALL( cudaMalloc( &cudaMemory, sizeof( float )*10 ), "malloc" );

	}
	void TearDown()
	{
		cudaFree( cudaMemory );
	}
};

TEST_F( ADeviceCommunicator, SetGetValueFloat )
{
	DeviceCommunicator<float>::setValue( cudaMemory, someIndex, someValue );

	ASSERT_FLOAT_EQ( someValue, DeviceCommunicator<float>::getValue( cudaMemory, someIndex ) );
}
