#include "gmock/gmock.h"
#include "LocalLink.h"
#include "su3/SU3Real12.h"
#include "su3/SU3Real18.h"
#include "su3/ParameterizationMediatorSU3_Real12_Real18.h"
#include "../cudatest/cuda_test_compare.h"
#include "../cudatest/cuda_test_outputtypes.h"

using namespace culgt;
using namespace ::testing;


namespace LocalLinkCu
{

struct DemonstrateCudaCompareTestFunc
{
	__host__ __device__ void operator()( CudaTestOutputFloat& out ) const
	{
#ifdef  __CUDA_ARCH__
		LocalLink<SU3Real18<float> > link18;
		link18.set(1, 2.42);
		out.result = link18.get( 1 );
#else
		LocalLink<SU3Real18<float> > link18;
		link18.set(1, 2.42);
		out.result = -1.;
#endif
	}
};

TEST( DemonstrateCudaCompareTest, ShouldBeFalseOnFunctionThatDoesDifferentThingsOnHostVsDevice )
{
	CudaTestOutputFloat out;
	cudaRunAndCompare( DemonstrateCudaCompareTestFunc(), out, false ); // skip last argument to check for equality
}

struct GetReturnsPreviouslySetValueFunc
{
	__host__ __device__ void operator()( CudaTestOutputFloat& out ) const
	{
		LocalLink<SU3Real18<float> > link18;
		link18.set(1, 2.42);
		out.result = link18.get( 1 );
#ifndef  __CUDA_ARCH__
		EXPECT_FLOAT_EQ( 2.42, out.result );
#endif
	}
};

TEST( ALocalLinkOnDevice, GetReturnsPreviouslySetValue )
{
	CudaTestOutputFloat out;
	cudaRunAndCompare( GetReturnsPreviouslySetValueFunc(), out );
}

struct OperatorAssignCopiesForSameParamTypeFunc
{
	__host__ __device__ void operator()( CudaTestOutputFloat& out ) const
	{
		LocalLink<SU3Real12<float> > link12_1;
		LocalLink<SU3Real12<float> > link12_2;
		link12_1.set(1, 2.42);
		link12_2 = link12_1;
		out.result = link12_2.get( 1 );
#ifndef  __CUDA_ARCH__
		EXPECT_FLOAT_EQ( 2.42, out.result );
#endif
	}
};

TEST( ALocalLinkOnDevice, OperatorAssignCopiesForSameParamType )
{
	CudaTestOutputFloat out;
	cudaRunAndCompare( OperatorAssignCopiesForSameParamTypeFunc(), out );
}

struct OperatorAssignCopiesFromReal18ToReal12Func
{
	__host__ __device__ void operator()( CudaTestOutputFloat& out ) const
	{
		LocalLink<SU3Real12<float> > link12;
		LocalLink<SU3Real18<float> > link18;

		link18.set(1, 2.42);
		link12 = link18;
		out.result = link12.get(1);
#ifndef  __CUDA_ARCH__
		EXPECT_FLOAT_EQ( 2.42, out.result );
#endif
	}
};

TEST( ALocalLinkOnDevice, OperatorAssignCopiesFromReal18ToReal12 )
{

	CudaTestOutputFloat out;
	cudaRunAndCompare( OperatorAssignCopiesFromReal18ToReal12Func(), out );
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

TEST( ALocalLinkOnDevice, OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLine )
{
	CudaTestOutputFloat out;
	cudaRunAndCompare( OperatorAssignCopiesFromReal12ToReal18AndReconstructsThirdLineFunc(), out );
}

}
