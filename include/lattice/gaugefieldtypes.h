#ifndef GAUGEFIELDTYPE_H
#define GAUGEFIELDTYPE_H

namespace culgt
{

namespace gaugefieldtype
{
	enum GaugefieldType {LINEAR, LOGARITHMIC};

	inline std::istream& operator>>(std::istream& in, gaugefieldtype::GaugefieldType& t)
	{
	    std::string token;
	    in >> token;
	    if ( boost::iequals(token, "LINEAR" ) )
	        t = LINEAR;
	    else if (boost::iequals(token, "LOGARITHMIC" ))
	        t = LOGARITHMIC;
	    return in;
	}

	inline std::ostream& operator<<(std::ostream& out, gaugefieldtype::GaugefieldType t)
	{
	    std::string token;
	    if (t == LINEAR)
	        token = "LINEAR";
	    else if (t == LOGARITHMIC)
	    	token = "LOGARITHMIC";
	    out << token;
	    return out;
	}
}

}

#endif
