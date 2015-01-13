#include <iostream>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/transform.hpp>

using namespace std;
namespace mpl = boost::mpl;

typedef mpl::vector_c< int, 4, 8 > aIntList;
typedef mpl::vector<mpl::string<'4'>, mpl::string<'8'> > aStringList;

typedef typename mpl::reverse_fold<aStringList,mpl::string<>,mpl::copy< mpl::copy<mpl::_1, mpl::back_inserter<mpl::string<'/'> > >, mpl::back_inserter<mpl::_2> > >::type constring;

template<typename theInt> struct TO_STRING
{
	typedef typename mpl::push_back< mpl::string<>, mpl::char_<'0'+theInt::value> >::type type;
//	typedef typename TO_STRING<typename theInt ::value>::type type;
};

template<int value> struct TO_STRING<mpl::int_<value> >
{
	typedef typename mpl::push_back< mpl::string<>, mpl::char_<'1'+value> >::type type;
};


int main()
{
	string test = mpl::c_str<constring>::value;
	cout << test << endl;


	// convert an int to a string
	typedef mpl::int_<4> my4;

	typedef mpl::transform< aIntList, TO_STRING<mpl::_> >::type convertedStringList;

	typedef typename mpl::reverse_fold<convertedStringList,mpl::string<>,mpl::copy< mpl::copy<mpl::_1, mpl::back_inserter<mpl::string<'/'> > >, mpl::back_inserter<mpl::_2> > >::type stringlist;
	cout << mpl::c_str<stringlist>::value << endl;

}
