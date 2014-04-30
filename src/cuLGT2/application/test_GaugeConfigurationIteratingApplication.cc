
#include "GaugeConfigurationIteratingApplication.h"
#include "gmock/gmock.h"
#include "../lattice/testhelper_pattern_stub.h"

using namespace testing;
using namespace culgt;

template<typename PatternType> class LinkFileStub: public LinkFile<PatternType>
{
public:
	LinkFileStub( LatticeDimension<4> dim ) : LinkFile<PatternType>::LinkFile( dim ) {};
	MOCK_METHOD0(loadImplementation,void());
	MOCK_METHOD0(saveImplementation,void());
};

class GaugeConfigurationIteratingApplicationMock: public GaugeConfigurationIteratingApplication<PatternStub<float>, LinkFileStub<PatternStub<float> > >
{
public:
	GaugeConfigurationIteratingApplicationMock( LatticeDimension<4> dim, FileIterator fileiterator ) : 									GaugeConfigurationIteratingApplication<PatternStub<float>, LinkFileStub<PatternStub<float> > >( dim, fileiterator, NULL ){}
	GaugeConfigurationIteratingApplicationMock( LatticeDimension<4> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : 	GaugeConfigurationIteratingApplication<PatternStub<float>, LinkFileStub<PatternStub<float> > >( dim, fileiterator, programOptions ){}
	MOCK_METHOD0(iterate,void());
};

TEST( AGaugeConfigurationIteratingApplication, RunCallsIterateNTimes )
{
	const int N = 10;
	LatticeDimension<4> dim( 4,4,4,4 );
	FileIterator fileiterator( "out_", ".dat", 4, 0, N-1, 1 );
	GaugeConfigurationIteratingApplicationMock app( dim, fileiterator );

	EXPECT_CALL( app, iterate() ).Times(N);

	app.run();
}

TEST( AGaugeConfigurationIteratingApplication, ConstructorWithCommandLineParametersAndBoostProgramOptions )
{
	const int argc_test = 7;
	const char* argv_test[argc_test] = {"NameOfProgram","--fbasename","test_","--fending",".test", "--nconf", "10"};

	GaugeConfigurationIteratingApplicationMock* app = GaugeConfigurationIteratingApplicationMock::init<GaugeConfigurationIteratingApplicationMock>( argc_test, argv_test );

	EXPECT_CALL( *app, iterate() ).Times(10);

	app->run();

	app->destroy();
}

class GaugeConfigurationIteratingApplicationWithFixedSettings: public Test
{
public:
	static const int N = 10;
	LatticeDimension<4> dim;
	FileIterator fileiterator;
	GaugeConfigurationIteratingApplicationMock app;

	GaugeConfigurationIteratingApplicationWithFixedSettings() : dim(4,4,4,4), fileiterator( "test_object_", ".dat", 4, 0, N-1, 1 ), app( dim, fileiterator )
	{
	}
};

TEST_F( GaugeConfigurationIteratingApplicationWithFixedSettings, LoadDelegatesLoadToLinkFileObject )
{
	EXPECT_CALL( app.getLinkFile(), loadImplementation() );
	app.load();
}

TEST_F( GaugeConfigurationIteratingApplicationWithFixedSettings, SaveDelegatesSaveToLinkFileObject )
{
	EXPECT_CALL( app.getLinkFile(), saveImplementation() );
	app.save( "test" );
}

