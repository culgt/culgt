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
typedef mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;

template<typename UseTexture, typename ThreadsPerSite> struct Printer2
{
	template<typename T> static void exec( T object )
	{
		cout << UseTexture::value << "/" << ThreadsPerSite::value << endl;
	}
};

struct UniqueId
{
};

int main()
{
	typedef RuntimeChooser<UniqueId,Printer2<_,_> > MyChooser;

	SequenceRunnerFrontend<MyChooser,useTextureSequence,threadsPerSiteSequence> test;

	for( vector<RuntimeChooserOption>::iterator it = MyChooser::begin(); it != MyChooser::end(); ++it )
	{
		cout << it->name << ": ";
		test.run( it->id );
	}
}
