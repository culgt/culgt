/*
 * filetypes.h
 *
 *  Created on: Jun 22, 2015
 *      Author: vogt
 */

#ifndef FILETYPES_H_
#define FILETYPES_H_
#include <boost/algorithm/string.hpp>

using boost::algorithm::iequals;

namespace culgt
{

enum LinkFileType {HEADERONLY, VOGT, HIREP, ILDG };

inline std::istream& operator>>(std::istream& in, LinkFileType& t)
{
    std::string token;
    in >> token;
    if ( boost::iequals(token, "HEADERONLY" ) )
        t = HEADERONLY;
    else if (boost::iequals(token, "VOGT" ))
        t = VOGT;
    else if (boost::iequals(token, "HIREP" ) )
    	t = HIREP;
    else if (boost::iequals(token, "ILDG" ) )
    	t = ILDG;
    return in;
}

inline std::ostream& operator<<(std::ostream& out, LinkFileType t)
{
    std::string token;
    if (t == HEADERONLY)
        token = "HEADERONLY";
    else if (t == VOGT)
    	token = "VOGT";
    else if (t == HIREP )
    	token = "HIREP";
    else if (t == ILDG )
    	token = "ILDG";
    out << token;
    return out;
}



}


#endif /* FILETYPES_H_ */
