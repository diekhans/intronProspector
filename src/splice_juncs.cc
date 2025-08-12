/*  splice_juncs.cc -- operations on junctions.

    Copyright (c) 2018-2025, Mark Diekhans, University of California, Santa Cruz
*/
#include "splice_juncs.hh"
#include <set>
#include <stdexcept>

// add all versions of one of  canonical junctions to set
static void add_canonical(const string& splice_sites,
                          int pos,
                          set<string>& canonicals) {
    canonicals.insert(splice_sites);
    string splice_sites_mod = splice_sites;
    for (int i = pos; i < splice_sites_mod.size(); i++) {
        splice_sites_mod[i] = ::toupper(splice_sites_mod[i]);
        add_canonical(splice_sites_mod, i + 1, canonicals);
        splice_sites_mod[i] = ::tolower(splice_sites_mod[i]);
        add_canonical(splice_sites_mod, i + 1, canonicals);
    }
}


// build all canonical junctions regardless of case
static void build_canonical(set<string>& canonicals) {
    add_canonical("GT/AG", 0, canonicals);
    add_canonical("GC/AG", 0, canonicals);
    add_canonical("AT/AC", 0, canonicals);
}

// is this canonical, junctions in the form GT/AG
bool is_canonical(const string& junctions) {
    static set<string> canonicals;  // has all combinations of upper and lower case
    if (canonicals.size() == 0) {
        build_canonical(canonicals);
    }
    return canonicals.count(junctions) > 0;
}

// parse JunctionFilter string
JunctionFilter junction_filter_parse(const string& spec) {
    if (spec == "all") {
        return ALL_SJ_FILTER;
    } else if (spec == "canon") {
        return CANON_SJ_FILTER;
    } else {
        throw invalid_argument("invalid splice junction filter '" + spec + "', expected 'all' or 'canon'");
    }
}

// check junction against filters, junctions in the form GT/AG
bool junction_filter_check(JunctionFilter junction_filter,
                           const string& splice_sites) {
    if (splice_sites.size() == 0) {
        return true;
    } else if (junction_filter == CANON_SJ_FILTER) {
        return is_canonical(splice_sites);
    } else {
        return true;
    }
}
