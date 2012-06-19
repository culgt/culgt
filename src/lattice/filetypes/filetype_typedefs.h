/*
 * filetype_typedefs.h
 *
 *  Created on: Jun 19, 2012
 *      Author: vogt
 */

#ifndef FILETYPE_TYPEDEFS_H_
#define FILETYPE_TYPEDEFS_H_

#include <boost/algorithm/string.hpp>

enum FileType {PLAIN, HEADERONLY, VOGT};

std::istream& operator>>(std::istream& in, FileType& t)
{
    std::string token;
    in >> token;
    if ( boost::iequals(token, "PLAIN" ) )
        t = PLAIN;
    else if (boost::iequals(token, "HEADERONLY" ))
        t = HEADERONLY;
    else if (boost::iequals(token, "VOGT" ) )
    	t = VOGT;
    return in;
}

std::ostream& operator<<(std::ostream& out, FileType& t)
{
    std::string token;
    if (t == PLAIN)
        token = "PLAIN";
    else if (t == HEADERONLY)
    	token = "HEADERONLY";
    else if (t == VOGT )
    	token = "VOGT";
    out << token;
    return out;
}


#endif /* FILETYPE_TYPEDEFS_H_ */
