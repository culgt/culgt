/*
 * filetypes.h
 *
 *  Created on: Jun 22, 2015
 *      Author: vogt
 */

#ifndef FILETYPES_H_
#define FILETYPES_H_

namespace culgt
{

namespace LinkFileType
{

enum FileType {DEFAULT,HEADERONLY, VOGT, HIREP, ILDG, NERSC, MDP};

inline std::istream& operator>>(std::istream& in, FileType& t)
{
    std::string token;
    in >> token;
    if ( token.compare( "HEADERONLY" ) == 0 )
        t = HEADERONLY;
    else if ( token.compare( "VOGT" ) == 0 )
        t = VOGT;
    else if ( token.compare( "HIREP" ) == 0 )
    	t = HIREP;
    else if ( token.compare( "ILDG" ) == 0 )
    	t = ILDG;
    else if ( token.compare( "NERSC" ) == 0 )
    	t = NERSC;
    else if ( token.compare( "MDP" ) == 0 )
    	t = MDP;
    else if ( token.compare( "DEFAULT" ) == 0 )
    	t = DEFAULT;
    return in;
}

inline std::ostream& operator<<(std::ostream& out, FileType t)
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
    else if (t == NERSC )
    	token = "NERSC";
    else if (t == DEFAULT )
    	token = "DEFAULT";
    out << token;
    return out;
}

}

}


#endif /* FILETYPES_H_ */
