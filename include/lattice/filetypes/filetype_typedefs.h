/**
 * filetype_typedefs.h
 *
 *  Created on: May 8, 2014
 *      Author: vogt
 */

#ifndef FILETYPE_TYPEDEFS_H_
#define FILETYPE_TYPEDEFS_H_

namespace culgt
{

enum ReinterpretReal {STANDARD, DOUBLE, FLOAT};

inline std::istream& operator>>(std::istream& in, ReinterpretReal& t)
{
    std::string token;
    in >> token;
    if ( token.compare( "STANDARD" ) == 0 )
        t = STANDARD;
    else if ( token.compare(  "DOUBLE" ) == 0 )
        t = DOUBLE;
    else if ( token.compare(  "FLOAT" ) == 0 )
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
