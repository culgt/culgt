#include "gmock/gmock.h"
#include "application/FileIterator.h"

using namespace testing;
using namespace culgt;
using namespace std;

TEST( AFileIterator, GetFilenameConcatenatesPartsOfDefaultConstructor )
{
	FileIterator iterator;

	ASSERT_EQ( "out_0000.dat", iterator.getFilename() );
}

TEST( AFileIterator, GetFilenameConcatenatesForCustomNameAndEnding )
{
	FileIterator iterator( "myBase_", ".end");

	ASSERT_EQ( "myBase_0000.end", iterator.getFilename() );
}

TEST( AFileIterator, GetFilenameConcatenatesWithNumberformat2 )
{
	FileIterator iterator( "myBase_", ".end", 2 );

	ASSERT_EQ( "myBase_00.end", iterator.getFilename() );
}

TEST( AFileIterator, GetFilenameConcatenatesWithFullConstructor )
{
	FileIterator iterator( "myBase_", ".end", 2, 43, 45, 1 );

	ASSERT_EQ( "myBase_43.end", iterator.getFilename() );
}

TEST( AFileIterator, HasElementIsFalseAfterNextWithStandardConstructor )
{
	FileIterator iterator;
	iterator.next();

	ASSERT_FALSE( iterator.hasElement() );
}

TEST( AFileIterator, HasNextIsTrueIfThereIsElement )
{
	const int someVal = 2;
	const int start = 0;
	const int end = 1;
	FileIterator iterator( "", "", someVal, start, end );

	ASSERT_TRUE( iterator.hasElement() );
	iterator.next();
	ASSERT_TRUE( iterator.hasElement() );
	iterator.next();
	ASSERT_FALSE( iterator.hasElement() );
}

TEST( AFileIterator, NextReturnsCurrentFileNameAndIncrements )
{
	const int start = 0;
	const int end = 1;
	FileIterator iterator( "out_", ".dat", 4, start, end );

	ASSERT_EQ( "out_0000.dat", iterator.next() );
	ASSERT_EQ( "out_0001.dat", iterator.next() );
}

TEST( AFileIterator, GetFileNameWithAppendix )
{
	FileIterator iterator;

	ASSERT_EQ( "out_appendix_0000.dat", iterator.getFilename( "appendix_" ) );
}

