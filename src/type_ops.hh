/*  type_ops.hh -- operations on types.

    Copyright (c) 2019, Mark Diekhans, University of California, Santa Cruz
**/
#ifndef type_ops_hh
#define type_ops_hh
#include <sstream>
#include <algorithm>
#include <stdint.h>
#include <stdexcept>

inline std::string int_to_string(int n) {
    std::stringstream s1;
    s1 << n;
    return s1.str();
}

inline uint32_t string_to_uint(const std::string& s) {
    char *end;
    uint32_t n = strtoul(s.c_str(), &end, 10);
    if ((end == s.c_str()) || (*end != '\0')) {
        throw std::runtime_error("invalid integer \"" + s + "\"");
    }
    return n;
}

inline uint64_t string_to_uint64(const std::string& s) {
    char *end;
    uint64_t n = strtoull(s.c_str(), &end, 10);
    if ((end == s.c_str()) || (*end != '\0')) {
        throw std::runtime_error("invalid integer \"" + s + "\"");
    }
    return n;
}

inline float string_to_float(const std::string& s) {
    char *end;
    float n = strtof(s.c_str(), &end);
    if ((end == s.c_str()) || (*end != '\0')) {
        throw std::runtime_error("invalid float \"" + s + "\"");
    }
    return n;
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

inline std::string mk_coords_str(const std::string &chrom, int start, int end) {
    return chrom + ":" + int_to_string(start) + "-" + int_to_string(end);
}

#endif
