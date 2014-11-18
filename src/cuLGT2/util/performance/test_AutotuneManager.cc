#include "gmock/gmock.h"
#include "AutotuneManager.h"

using namespace testing;
using namespace culgt;
using namespace std;

class AAutotuneManager : public Test
{
public:
	AutotuneManager autotune;
};

TEST_F( AAutotuneManager, AddSingleAttributeGetIdentifier )
{
	string testAttribute = "someTestAttribute";

	autotune.addAttribute( testAttribute );

	ASSERT_EQ( "tune_" + testAttribute, autotune.getIdentifier() );
}

TEST_F( AAutotuneManager, TryLoadOptimalIdFromNonExistentFile )
{
	string testAttribute = "someNonExistingId";

	autotune.addAttribute( testAttribute );

	ASSERT_THROW( autotune.getOptimalId(), AutotuneManagerOptionNotAvailable );
}

/**
 * Ensure the file "tune_optimalid" exists and has number zero in it.
 */
TEST_F( AAutotuneManager, LoadsOptimalIdFromFile )
{
	string testAttribute = "optimalid";

	autotune.addAttribute( testAttribute );

	ASSERT_EQ( 0, autotune.getOptimalId() );
}

TEST_F( AAutotuneManager, WritesOptimalIdToFile )
{
	string testAttribute = "writtenOptimalid";
	int optimalId = time( NULL );

	autotune.addAttribute( testAttribute );
	autotune.writeOptimalId( optimalId );

	ASSERT_EQ( optimalId, autotune.getOptimalId() );
}

