#include <iostream>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/placeholders.hpp>

#include "../../cudacommon/LaunchBounds.h"
#include "SequenceRunner.h"
#include "RuntimeChooser.h"

namespace mpl = boost::mpl;
using mpl::placeholders::_;

typedef mpl::vector_c< int, 0, 1 > useTextureSequence;

template<typename UseTexture> struct Printer
{
	template<typename T> static void exec( T object )
	{
		cout << UseTexture::value << endl;
	}
};

struct UniqueId
{
};

int main()
{
	typedef RuntimeChooser<UniqueId,Printer<_> > MyChooser;

	SequenceRunnerFrontend<MyChooser,useTextureSequence> test;

	for( vector<RuntimeChooserOption>::iterator it = MyChooser::begin(); it != MyChooser::end(); ++it )
	{
		cout << it->name << ": ";
		test.run( it->id );
	}
}
