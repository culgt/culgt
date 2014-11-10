/*
 * RuntimeChooser.h
 *
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef RUNTIMECHOOSER_H_
#define RUNTIMECHOOSER_H_
#include <vector>
#include <boost/functional/hash.hpp>


using std::vector;

template<typename GiveAutoTuneClassHere, typename Run> class RuntimeChooser
{
public:
	static size_t id;
	static vector<size_t> ids;
	static bool init;

	static vector<size_t>::iterator begin()
	{
		return ids.begin();
	}

	static vector<size_t>::iterator end()
	{
		return ids.end();
	}

	template<typename T> void operator()(T) const // T should be a Sequence
	{
		boost::hash<std::string> string_hash;
		size_t hashed = string_hash(typeid(T).name());

//#ifdef __INTEL_COMPILER
		// seems to not support c++11 function typeid().hash_code()
//		boost::hash<std::string> string_hash;
//		size_t hashed = string_hash(typeid(T).name());
//#else
//		size_t hashed = typeid(T).hash_code(); // the standard does not demand that this does not change between runs... (maybe same for boost::hash())
//#endif

		if( init )
		{
			ids.push_back( hashed );
		}
		else if( id ==  hashed )
		{
			typedef mpl::unpack_args<typename mpl::lambda<Run>::type > g;
			typename mpl::apply< g, typename T::type >::type t;
			t();
		}
	}
};

template<typename T0, typename T1> size_t RuntimeChooser<T0,T1>::id;
template<typename T0, typename T1> vector<size_t> RuntimeChooser<T0,T1>::ids;
template<typename T0, typename T1> bool RuntimeChooser<T0,T1>::init = false;


#endif /* RUNTIMECHOOSER_H_ */
