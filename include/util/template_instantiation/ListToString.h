/*
 * ListToString.h
 *
 *  Created on: Jan 13, 2015
 *      Author: vogt
 */

#ifndef LISTTOSTRING_H_
#define LISTTOSTRING_H_

#include <sstream>
#include <string>
#include <boost/mpl/for_each.hpp>

namespace culgt
{
namespace mpl
{
	template<typename UniqueID> struct StringContainer
	{
		static std::stringstream data;
	};
	template<typename UniqueID> std::stringstream StringContainer<UniqueID>::data;

	template<typename UniqueID> struct to_string
	{
		template<int T> void operator()(boost::mpl::int_<T>) const
		{
			StringContainer<UniqueID>::data << T << "/";
		}
		template<int T> void operator()(boost::mpl::integral_c<int, T>) const
		{
			StringContainer<UniqueID>::data << T << "/";
		}
		template<typename T> void operator()(T) const
		{
			StringContainer<UniqueID>::data << T::printable() << "/";
		}
	};

	template<typename List> class ListToString
	{
	public:
		static std::string getString()
		{
			boost::mpl::for_each<List>(to_string<List>());
			std::string str = StringContainer<List>::data.str();
			return str.substr(0,str.size()-1);
		}
	};
}
}

#endif /* LISTTOSTRING_H_ */
