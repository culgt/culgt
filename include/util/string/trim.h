/*
 * Using boost's trim gives an error in cuda <= 7.0 due to a bug
 */

#ifndef TRIM_H_
#define TRIM_H_

#include <string>

namespace culgt
{

// trim from start
static inline std::string& ltrim(std::string &s)
{
	size_t start = s.find_first_not_of( " \t\n\v\f\r" );
	s.erase(0, start );
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
		size_t end = s.find_last_not_of( " \t\n\v\f\r" );
        s.erase(end+1, s.size()-end);
        return s;
}
// trim from both ends
static inline std::string& trim(std::string &s) {
        return ltrim(rtrim(s));
}

}



#endif /* TRIM_H_ */
