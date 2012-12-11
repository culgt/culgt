/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 */

#ifndef FILETYPE_TYPEDEFS_H_
#define FILETYPE_TYPEDEFS_H_

#include <boost/algorithm/string.hpp>


enum ReinterpretReal {STANDARD, DOUBLE, FLOAT};

std::istream& operator>>(std::istream& in, ReinterpretReal& t)
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

std::ostream& operator<<(std::ostream& out, ReinterpretReal& t)
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
