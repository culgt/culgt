#include "gmock/gmock.h"
#include "util/string/trim.h"

using namespace testing;
using namespace culgt;
using namespace std;

TEST( StringTrim, Left )
{
	string str = " \t trim\t ";

	ltrim(str);

	ASSERT_STREQ( "trim\t ", str.c_str());
}

TEST( StringTrim, Right )
{
	string str = " \t trim\t ";

	rtrim(str);

	ASSERT_STREQ( " \t trim", str.c_str());
}

TEST( StringTrim, Both )
{
	string str = " \t trim\t ";

	trim(str);

	ASSERT_STREQ( "trim", str.c_str());
}
