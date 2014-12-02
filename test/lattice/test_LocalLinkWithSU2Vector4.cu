#include "gmock/gmock.h"
#include "lattice/LocalLink.h"
#include "lattice/parameterization_types/SU2Vector4.h"
#include "lattice/parameterization_types/SUNRealFull.h"
#include "lattice/parameterization_types/ParameterizationMediatorSU2_Vector4_Real8.h"

using namespace culgt;
using namespace ::testing;




class ALocalLinkWithSU2Vector4: public Test
{
public:
	LocalLink<SUNRealFull<2,float> > linkReal8;
	LocalLink<SU2Vector4<float> > linkVector4;
	LocalLink<SU2Vector4<float> > linkVector4_2;

	float4 aFloat4;
	void SetUp()
	{
		aFloat4.x = 1.1;
		aFloat4.y = 2.2;
		aFloat4.z = 3.3;
		aFloat4.w = 4.4;
		linkReal8.zero();
		linkVector4.zero();
		linkVector4_2.zero();
	}
};

TEST_F( ALocalLinkWithSU2Vector4, GetReturnsPreviouslySetValue )
{
	linkVector4.set( 0, aFloat4 );
	ASSERT_FLOAT_EQ( aFloat4.y, linkVector4.get(0).y );
}

TEST_F( ALocalLinkWithSU2Vector4, OperatorAssignCopiesForSameParamType )
{
	linkVector4.set( 0, aFloat4 );

	linkVector4_2 = linkVector4;

	ASSERT_FLOAT_EQ( aFloat4.y, linkVector4_2.get(0).y );
}

TEST_F( ALocalLinkWithSU2Vector4, OperatorAssignCopiesFromVector4ToRealFull )
{
	linkVector4.set( 0, aFloat4 );

	linkReal8 = linkVector4;

	ASSERT_FLOAT_EQ( -aFloat4.y, linkReal8.get(7) );
}

TEST_F( ALocalLinkWithSU2Vector4, OperatorAssignCopiesFromRealFullToVector4 )
{
	linkReal8.set( 1, 1.5 );

	linkVector4 = linkReal8;

	ASSERT_FLOAT_EQ( 1.5, linkVector4.get(0).y );
}
