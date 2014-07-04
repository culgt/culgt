/**
 * filetype_typedefs.h
 *
 *  Created on: May 8, 2014
 *      Author: vogt
 */

#ifndef FILETYPE_TYPEDEFS_H_
#define FILETYPE_TYPEDEFS_H_

#include <boost/algorithm/string.hpp>

namespace culgt
{

enum ReinterpretReal {STANDARD, DOUBLE, FLOAT};

inline std::istream& operator>>(std::istream& in, ReinterpretReal& t)
{
    std::string token;
    in >> token;
    if ( boost::iequals(token, "STANDARD" ) )
        t = STANDARD;
    else if (boost::iequals(token, "DOUBLE" ))
        t = DOUBLE;
    else if (boost::iequals(token, "FLOAT" ) )
    	t = FLOAT;
    return in;
}

inline std::ostream& operator<<(std::ostream& out, ReinterpretReal t)
{
    std::string token;
    if (t == STANDARD)
        token = "STANDARD";
    else if (t == DOUBLE)
    	token = "DOUBLE";
    else if (t == FLOAT )
    	token = "FLOAT";
    out << token;
    return out;
}

}


#endif /* FILETYPE_TYPEDEFS_H_ */
