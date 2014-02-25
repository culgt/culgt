#include "gmock/gmock.h"
#include "LocalLink.h"
#include "su3/SU3Real12.h"
#include "su3/SU3Real18.h"
#include "su3/ParameterizationMediatorSU3_Real12_Real18.h"
#include "../cudatest/cuda_test_compare.h"
#include "../cudatest/cuda_gtest_plugin.h"
#include "../cudatest/cuda_test_outputtypes.h"

using namespace culgt;
using namespace ::testing;


namespace LocalLinkCu
{

CUDA_TEST( ALocalLinkOnDevice, GetReturnsPreviouslySetValue )
{
	LocalLink<SU3Real18<float> > link18;
	link18.set(1, 2.42);
	ASSERT_FLOAT_EQ( 2.42f, link18.get( 1 ) );
}

CUDA_TEST( ALocalLinkOnDevice, OperatorAssignCopiesForSameParamType )
{
	float someValue = 2.42;
	LocalLink<SU3Real12<float> > link12_1;
	LocalLink<SU3Real12<float> > link12_2;
	link12_1.set(1, someValue);

	link12_2 = link12_1;

	ASSERT_FLOAT_EQ( someValue, link12_2.get( 1 ) );
}

CUDA_TEST( ALocalLinkOnDevice, OperatorAssignCopiesFromReal18ToReal12 )
{
	float someValue = 2.42;
	LocalLink<SU3Real12<float> > link12;
	LocalLink<SU3Real18<float> > link18;

	link18.set(1, someValue );
	link12 = link18;
	ASSERT_FLOAT_EQ( someValue, link12.get( 1 ) );
}

struct OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLineFunc
{
	__host__ __device__ void operator()( CudaTestOutputFloat& out ) const
	{
		LocalLink<SU3Real12<float> > link12;
		LocalLink<SU3Real18<float> > link18;

		link12.zero();
		link12.set(0, 1.);
		link12.set(8, 1.);

		link18 = link12;
		out.result = link18.get(16);
#ifndef  __CUDA_ARCH__
		EXPECT_FLOAT_EQ( 1., out.result );
#endif
	}
};

CUDA_TEST( ALocalLinkOnDevice, OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLine )
{
	LocalLink<SU3Real12<float> > link12;
	LocalLink<SU3Real18<float> > link18;
	link12.zero();
	link12.set(0, 1.);
	link12.set(8, 1.);

	link18 = link12;

	ASSERT_FLOAT_EQ( 1.f, link18.get( 16 ) );
}

}
