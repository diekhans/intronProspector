/*  type_ops.hh -- operations on types.

    Copyright (c) 2019, Mark Diekhans, University of California, Santa Cruz
**/
#ifndef type_ops_hh
#define type_ops_hh
#include <sstream>
#include <algorithm>

// convert an integer to a string (newer version of C++ have std::to_string)
inline std::string int_to_string(int num) {
    std::stringstream s1;
    s1 << num;
    return s1.str();
}

inline std::string to_upper(const std::string& str) {
    std::string ustr(str);
    std::transform(ustr.begin(), ustr.end(), ustr.begin(), ::toupper);
    return ustr;
}

inline std::string to_lower(const std::string& str) {
    std::string lstr(str);
    std::transform(lstr.begin(), lstr.end(), lstr.begin(), ::tolower);
    return lstr;
}

#endif
