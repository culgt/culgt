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
#include <typeinfo>
#include <boost/mpl/at.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_float.hpp>
#include "ListToString.h"
#include "RuntimeChooserOption.h"

using std::vector;

template<typename AutoTuneClass, typename Run> class RuntimeChooser
{
public:
	static size_t id;
	static vector<RuntimeChooserOption> options;
	static bool init;
	static AutoTuneClass* object;
	typedef void (*FPTR)(AutoTuneClass*);
	static FPTR run;
	static bool functionPtrIsSet;

	static vector<RuntimeChooserOption>::iterator begin()
	{
		return options.begin();
	}

	static vector<RuntimeChooserOption>::iterator end()
	{
		return options.end();
	}

	template<typename T> void operator()(T) const // T should be a Sequence
	{
		boost::hash<std::string> string_hash;
		size_t hashed = string_hash(typeid(T).name());

		if( init )
		{
			RuntimeChooserOption option;
			option.id = hashed;
			option.name = culgt::mplextension::ListToString<T>::getString();

			options.push_back( option );
		}
		else if( id == hashed )
		{
			typedef mpl::unpack_args<typename mpl::lambda<Run>::type > g;
			run = mpl::apply< g, typename T::type >::type::exec;
		}
	}
};

template<typename T0, typename T1> size_t RuntimeChooser<T0,T1>::id;
template<typename T0, typename T1> vector<RuntimeChooserOption> RuntimeChooser<T0,T1>::options;
template<typename T0, typename T1> bool RuntimeChooser<T0,T1>::init = false;
template<typename T0, typename T1> T0* RuntimeChooser<T0,T1>::object;
template<typename T0, typename T1> typename RuntimeChooser<T0,T1>::FPTR RuntimeChooser<T0,T1>::run;
template<typename T0, typename T1> bool RuntimeChooser<T0,T1>::functionPtrIsSet = false;


#endif /* RUNTIMECHOOSER_H_ */
