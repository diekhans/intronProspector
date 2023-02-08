/*
 * string to number conversion functions
 */
#ifndef STRCONVERT_H
#define STRCONVERT_H
#include <string>
#include <exception>
using namespace std;

/*
 * Convert a string to a float.
 */
inline float toFloat(const char* str) {
    char *endPtr;
    errno = 0;
    float num = strtof(str, &endPtr);
    if ((endPtr == str) || (*endPtr != '\0') || (errno != 0)) {
        if (errno != 0) {
            throw runtime_error("Float out of range \"" + string(str) + "\"");
        } else {
            throw runtime_error("Invalid float \"" + string(str) + "\"");
        }
    }
    return num;
}

/*
 * Convert a string to an unsigned.
 */
inline unsigned toUnsigned(const char* str,
                           int base = 10) {
    char *endPtr;
    errno = 0;
    unsigned long lnum = strtoul(str, &endPtr, base);
    if ((endPtr == str) || (*endPtr != '\0')) {
        throw runtime_error("Invalid unsigned integer \"" + string(str) + "\"");
    }
     
    unsigned num = (unsigned)lnum;
    if ((errno != 0) || ((unsigned long)num != lnum)) {
        throw runtime_error("Unsigned integer out of range \"" + string(str) + "\"");
    }
    return num;
}

/*
 * Convert a string to a long long.
 */
inline long long toLongLong(const char* str) {
    char *endPtr;
    errno = 0;
    long long num = strtoll(str, &endPtr, 0);
    if ((endPtr == str) || (*endPtr != '\0')) {
        throw runtime_error("Invalid long long \"" + string(str) + "\"");
    }
    return num;
}


#endif
