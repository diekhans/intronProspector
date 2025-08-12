/*  splice_juncs.hh -- operations on junctions.

    Copyright (c) 2018-2025, Mark Diekhans, University of California, Santa Cruz

*/
#ifndef SPLICE_JUNCS_HH
#define SPLICE_JUNCS_HH
using namespace std;
#include <string>

// junction filtering categories
typedef enum {
    NULL_SJ_FILTER = 0,   // no filtering specified
    ALL_SJ_FILTER = 1,    // all junctions
    CANON_SJ_FILTER = 2,  // canonical junctions
} JunctionFilter;

// is this canonical, junctions in the form GT/AG
bool is_canonical(const string& junctions);

// parse JunctionFilter string
JunctionFilter junction_filter_parse(const string& spec);

// check junction against filters, junctions in the form GT/AG
bool junction_filter_check(JunctionFilter junction_filter,
                           const string& splice_sites);
#endif
