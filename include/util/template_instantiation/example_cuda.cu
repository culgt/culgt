/**
 * This is a technique to choose from a template instantiation at RUNTIME!
 *
 * Needs:
 *  a) -std=c++11
 *  b) boost::mpl
 *
 * Example: In our CUDA autotuner we need to try the kernel with a set of different template parameters and then choose the optimal one.
 * We have 3 different parameters that can take several values and we want to try all combinations.
 * 1) Make a list for each parameter with the values to try.
 * 2) Define a method to call and hide it in a templated struct. The template arguments are the parameters. The method is defined in the operator().
 * 3) Define a struct or whatever you want that can serve as a unique identifier of your RuntimeChooser.
 * 4) typedef the unique RuntimeChooser for convenience. Template arguments are (a) the unique identifier from (3) and (b) the struct that holds the method to call from (2)
 * 	  where the template arguments are replaced by an underscore "_".
 * 5) Instantiate the frontend. Pass the RuntimeChooser and all sequences in the order in which the parameters are used in the function structure (2).
 * 6) From the RuntimeChooser you can get an iterator with all possible choices as a hash number.
 * 7) The frontend instance run() can now be called with one of the possible hash numbers.
 */

#include <iostream>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/placeholders.hpp>

#include "../../cudacommon/LaunchBounds.h"
#include "SequenceRunner.h"
#include "RuntimeChooser.h"
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>

using mpl_::int_;
using mpl_::bool_;

//namespace mpl = boost::mpl;
//using mpl::placeholders::_;

/*
 * 1) make the lists
 */
typedef mpl::vector< LaunchBounds<1,2>, LaunchBounds<3,4> > launchBoundsSequence;
typedef mpl::vector_c< int, 4, 8 > threadsPerSiteSequence;
typedef mpl::vector_c< int, 0, 1 > useTextureSequence;


template<typename LaunchBounds, typename ThreadsPerSite, typename UseTexture> __global__ void someKernel( int* a )
{
	printf( "device: (%d,%d), %d, %d\n", LaunchBounds::maxThreadsPerBlock, LaunchBounds::minBlocksPerMultiprocessor, ThreadsPerSite::value, UseTexture::value );
	a[threadIdx.x] = threadIdx.x;
}

/*
 * 2) define the method to call: here we just print out the values.
 */
template<typename LaunchBounds, typename ThreadsPerSite, typename UseTexture> struct Printer
{
	template<typename T> static void exec( T* object )
	{
		someKernel<LaunchBounds,ThreadsPerSite,UseTexture><<<1,1>>>( object->devPtr );
		cout << "(" << LaunchBounds::maxThreadsPerBlock << "," << LaunchBounds::minBlocksPerMultiprocessor << "), " << ThreadsPerSite::value << ", " << UseTexture::value << endl;
	}
};

/*
 * 3) Define a unique identifier
 *    Background: The possible runtime choices in the RuntimeChooser are static member variables, if you have several RuntimeChooser in your program we need to distinguish them...
 */
struct UniqueId
{
	int* devPtr;

	void allocate()
	{
		cudaMalloc( &devPtr, sizeof(int) );
	}
};

int main()
{
	UniqueId obj;
	obj.allocate();

	// 4) typedef the RuntimeChooser, here we want to call the Printer template with three arguments: Printer<_,_,_>
	typedef RuntimeChooser<UniqueId,Printer<_,_,_> > MyChooser;

	MyChooser::object = &obj;

//	Printer<LaunchBounds<4,3>,int_<3>, bool_<true> >::exec( MyChooser::object );



//	 5) instantiate the frontend
	SequenceRunnerFrontend<MyChooser,launchBoundsSequence,threadsPerSiteSequence,useTextureSequence> test;

	// 6) get an iterator and loop over all possible choices.
	for( vector<RuntimeChooserOption>::iterator it = MyChooser::begin(); it != MyChooser::end(); ++it )
	{
		cout << it->name << ": ";
		test.run( it->id );
	}

	int hostPtr[1];
	cudaMemcpy( hostPtr, obj.devPtr, sizeof(int), cudaMemcpyDeviceToHost );

}
