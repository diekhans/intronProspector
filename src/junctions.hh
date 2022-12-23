/*  junctions_extractor.hh -- splice junction extraction.

    Copyright (c) 2018, Mark Diekhans, University of California, Santa Cruz

This code is derived from the regtools package available at:

  https://github.com/griffithlab/regtools

We thank them for providing this software under a flexible license:

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef JUNCTIONS_H
#define JUNCTIONS_H

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include "htslib/sam.h"

using namespace std;

// Categories of reads
typedef enum {
    SINGLE_MAPPED_READ = 0,  // confident single-mapped
    MULTI_MAPPED_READ = 1,   // confident multi-mapped
    UNSURE_READ = 2,         // discordant, mate-unmapped, NH flag not set
} ReadCategory;
static const unsigned READ_CATEGORY_MAX = UNSURE_READ;



// class used to map Junction objects
class JunctionKey {
    public:
    const string chrom;
    const uint32_t intron_start;
    const uint32_t intron_end;

    JunctionKey(const string& chrom,
                uint32_t intron_start,
                uint32_t intron_end):
        chrom(chrom), intron_start(intron_start), intron_end(intron_end) {
    }

    string to_string() const {
        stringstream s1;
        s1 << chrom << ":" << intron_start << "-" << intron_end;
        return s1.str();
    }
};

// comparison key
static inline bool operator<(const JunctionKey &jk1, const JunctionKey &jk2) {
    if (jk1.chrom < jk2.chrom) {
        return true;
    }
    if (jk1.chrom > jk2.chrom) {
        return false;
    }
    // Same chromosome
    if (jk1.intron_start < jk2.intron_start) {
        return true;
    }
    if (jk1.intron_start > jk2.intron_start) {
        return false;
    }
    // Same start
    if (jk1.intron_end < jk2.intron_end) {
        return true;
    }
    if (jk1.intron_end > jk2.intron_end) {
        return false;
    }
    // Same end
    return false;
}

// Data save for an intron
class Junction {
public:
    static float NULL_CONFIDENCE;
    
    const string chrom;
    uint32_t intron_start;
    uint32_t intron_end;
    uint32_t anchor_start;  // This is the intron_start - max overhang
    uint32_t anchor_end;    // This is the end + max overhang
    char strand;
    string splice_sites;    // splite sites in the form GT/AG if genome available
    uint32_t ijunc;         // used for generate BED names

    private:
    // Number of reads supporting the junction, by category
    uint64_t read_counts[READ_CATEGORY_MAX + 1];
    // For each read, the distance from the start of the read to the start of the intron.
    // Uses in calculating Shannon-Wiener Diversity Index
    std::vector<uint16_t> read_offsets;
    float confidence;   // computed in a lazy manner

    void lazy_get_confidence();
    float calculate_confidence();
    
    public:
    Junction(const string& chrom1, uint32_t intron_start1, uint32_t intron_end1,
             uint32_t anchor_start1, uint32_t anchor_end1,
             char strand1, const string& splice_sites, uint32_t ijunc):
        chrom(chrom1), intron_start(intron_start1), intron_end(intron_end1),
        anchor_start(anchor_start1), anchor_end(anchor_end1),
        strand(strand1), splice_sites(splice_sites),
        ijunc(ijunc), confidence(NULL_CONFIDENCE) {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = 0;
        }
    }
    Junction(const Junction& src):
        chrom(src.chrom), intron_start(src.intron_start), intron_end(src.intron_end),
        anchor_start(src.anchor_start), anchor_end(src.anchor_end),
        strand(src.strand), splice_sites(src.splice_sites), ijunc(src.ijunc),
        read_offsets(src.read_offsets), confidence(src.confidence) {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = src.read_counts[i];
        }
    }

    string get_description() const;
    
    void set_read_counts(ReadCategory read_category,
                         uint32_t count) {
        read_counts[read_category] = count;
    }

    void set_confidence(float confidence1) {
        confidence = confidence1;
    }

    void merge(const Junction& j1);
    
    // add a read to the counts and offsets
    void count_read(ReadCategory read_category,
                    uint32_t read_offset) {
        read_counts[read_category] += 1;
        read_offsets.push_back(read_offset);
    }

    // sum different read counts
    uint64_t total_read_count() const {
        uint64_t cnt = 0;
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            cnt += read_counts[i];
        }
        return cnt;
    }

    // get the confidence, computing for the first time
    float get_confidence() const {
        if (isnan(confidence)) {
            const_cast<Junction*>(this)->lazy_get_confidence();
        }
        return confidence;
    }

    // is this a canonical intron?
    bool is_canonical() const {
        return (splice_sites.size() > 0) and (splice_sites[0] >= 'A')
            and (splice_sites[0] <= 'Z');
    }
    
    // Print BED with anchors as blocks and intron as gap.
    void print_anchor_bed(ostream& out) const;

    // Print BED with intron as block 
    void print_intron_bed(ostream& out) const;

    // Print row to junction call TSV
    void print_junction_call_row(ostream& out) const;
};

// Vector of pointers to junctions
class JunctionVector: public vector<Junction*> {
    private:
    static bool junctions_lt(const Junction *j1,
                             const Junction *j2);
    static bool introns_lt(const Junction *j1,
                           const Junction *j2);

    public:
    // sort by introns start
    void sort_by_introns() {
        sort(begin(), end(), introns_lt);
    }

    // sort by anchor start
    void sort_by_anchors() {
        sort(begin(), end(), junctions_lt);
    }
};

// Table of junctions, indexed by key
class JunctionTable: public map<JunctionKey, Junction*> {
    public:
    ~JunctionTable() {
        clear();
    }

    // get all junctions
    JunctionVector get_junctions() {
        JunctionVector juncs;
        for (JunctionTable::iterator it = begin(); it != end(); it++) {
            juncs.push_back(it->second);
        }
        return juncs;
    }

    void clear() {
        for (JunctionTable::iterator it = begin(); it != end(); it++) {
            delete it->second;
        }
        map<JunctionKey, Junction*>::clear();
    }
};

// Print BED with anchors as blocks and intron as gap.
void print_anchor_bed(const JunctionVector& juncs,
                      float min_confidence_score,
                      ostream& out);

// Print BED with intron as block
void print_intron_bed(const JunctionVector& juncs,
                      float min_confidence_score,
                      ostream& out);

// Print header for junction call TSV
void print_junction_call_header(ostream& out);

// Print TSV with intron information
void print_intron_call_tsv(const JunctionVector& juncs,
                           float min_confidence_score,
                           ostream& out);

#endif


