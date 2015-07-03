/*
 * culgt_compat.h
 *
 *  Created on: Jul 1, 2015
 *      Author: vogt
 */

#ifndef CULGT_COMPAT_H_
#define CULGT_COMPAT_H_

#if __cplusplus >= 201103L
#define CULGT_OVERRIDE override
#define CULGT_FINAL final
#else
#define CULGT_OVERRIDE
#define CULGT_FINAL
#endif

#include <string>
#include <sstream>
namespace culgt
{
	inline std::string to_string( int val )
	{
#if __cplusplus >= 201103L
		return std::to_string( val );
#else
		std::ostringstream oss;
		oss << val;
		return oss.str();
#endif
	}
}


#endif /* CULGT_COMPAT_H_ */
