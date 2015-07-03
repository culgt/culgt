
#include "application/GaugeConfigurationIteratingApplication.h"
#include "gmock/gmock.h"
#include "../lattice/testhelper_pattern_stub.h"

using namespace testing;
using namespace culgt;

template<typename PatternType> class LinkFileStub: public LinkFile<PatternType>
{
public:
	LinkFileStub( LatticeDimension<4> dim, ReinterpretReal reinterpret ) : LinkFile<PatternType>::LinkFile( dim, reinterpret ) {};
	MOCK_METHOD0(loadImplementation,void());
	MOCK_METHOD0(saveImplementation,void());
};

class GaugeConfigurationIteratingApplicationMock: public GaugeConfigurationIteratingApplication<PatternStub<float> >
{
public:
	GaugeConfigurationIteratingApplicationMock( LatticeDimension<4> dim, FileIterator fileiterator ) : 									GaugeConfigurationIteratingApplication<PatternStub<float> >( dim, fileiterator, NULL ){}
	GaugeConfigurationIteratingApplicationMock( LatticeDimension<4> dim, FileIterator fileiterator, ProgramOptions* programOptions ) : 	GaugeConfigurationIteratingApplication<PatternStub<float> >( dim, fileiterator, programOptions ){}
	MOCK_METHOD0(iterate,void());
	MOCK_METHOD0(setup,void());
	MOCK_METHOD0(teardown,void());
};

TEST( AGaugeConfigurationIteratingApplication, RunCallsSetup_IterateNTimes_Teardown )
{
	const int N = 10;
	LatticeDimension<4> dim( 4,4,4,4 );
	FileIterator fileiterator( "out_", ".dat", 4, 0, N-1, 1 );
	GaugeConfigurationIteratingApplicationMock app( dim, fileiterator );


	EXPECT_CALL( app, setup() );
	EXPECT_CALL( app, iterate() ).Times(N);
	EXPECT_CALL( app, teardown() );

	app.run();
}

TEST( AGaugeConfigurationIteratingApplication, ConstructorWithCommandLineParametersAndBoostProgramOptions )
{
	const int argc_test = 7;
	const char* argv_test[argc_test] = {"NameOfProgram","--fbasename","test_","--fextension",".test", "--nconf", "10"};

	GaugeConfigurationIteratingApplicationMock* app = GaugeConfigurationIteratingApplicationMock::init<GaugeConfigurationIteratingApplicationMock>( argc_test, argv_test );

	EXPECT_CALL( *app, setup() );
	EXPECT_CALL( *app, iterate() ).Times(10);
	EXPECT_CALL( *app, teardown() );

	app->run();

	app->destroy();
}

//class GaugeConfigurationIteratingApplicationWithFixedSettings: public Test
//{
//public:
//	static const int N = 10;
//	LatticeDimension<4> dim;
//	FileIterator fileiterator;
////	GaugeConfigurationIteratingApplicationMock app;
//
//	static const int argc_test = 1;
//	const char* argv_test[argc_test];
//
//	GaugeConfigurationIteratingApplicationMock* app;
//
//	GaugeConfigurationIteratingApplicationWithFixedSettings() : dim(4,4,4,4), fileiterator( "test_object_", ".dat", 4, 0, N-1, 1 )//, app( dim, fileiterator )
//	{
//	}
//
//	~GaugeConfigurationIteratingApplicationWithFixedSettings()
//	{
////		app->destroy();
//	}
//
//	void TearDown()
//	{
//	}
//};
//
//TEST_F( GaugeConfigurationIteratingApplicationWithFixedSettings, LoadDelegatesLoadToLinkFileObject )
//{
//	app = GaugeConfigurationIteratingApplicationMock::init<GaugeConfigurationIteratingApplicationMock>( argc_test, argv_test );
//
//	EXPECT_CALL( app->getLinkFile(), loadImplementation() );
//	app->load();
//
////	Mock::VerifyAndClearExpectations(&app->getLinkFile());
//}

//TEST_F( GaugeConfigurationIteratingApplicationWithFixedSettings, SaveDelegatesSaveToLinkFileObject )
//{
//	EXPECT_CALL( app.getLinkFile(), saveImplementation() );
//	app.save( "test" );
//}

