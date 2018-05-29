/*  junctions_extractor.h -- splice junction extraction.

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

#ifndef JUNCTIONS_EXTRACTOR_H
#define JUNCTIONS_EXTRACTOR_H

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include "htslib/sam.h"

using namespace std;

// class used to map Junction objects
class JunctionKey {
    public:
    const string& chrom;
    const uint32_t intron_start;
    const uint32_t intron_end;
    const char strand;

    JunctionKey(const string& chrom,
                uint32_t intron_start,
                uint32_t intron_end,
                char strand):
        chrom(chrom), intron_start(intron_start), intron_end(intron_end), strand(strand) {
    }
};

// comparison key
static inline bool operator<(const JunctionKey &jk1, const JunctionKey &jk2) {
    if (jk1.chrom < jk2.chrom){
        return true;
    }
    if (jk1.chrom > jk2.chrom){
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
    return jk1.strand < jk2.strand;
}

// Data save for an intron
class Junction {
public:
    const string& chrom;
    uint32_t intron_start;
    uint32_t intron_end;
    char strand;
    // Number of reads supporting the junction
    unsigned int read_count;
    // This is the intron_start - max overhang
    uint32_t anchor_start;
    // This is the end + max overhang
    uint32_t anchor_end;

    Junction(const string& chrom1, uint32_t intron_start1, uint32_t intron_end1,
             uint32_t anchor_start1, uint32_t anchor_end1,
             char strand1):
        chrom(chrom1), intron_start(intron_start1), intron_end(intron_end1),
        anchor_start(anchor_start1), anchor_end(anchor_end1),
        strand(strand1), read_count(0) {
    }
    Junction(const Junction& src):
        chrom(src.chrom), intron_start(src.intron_start), intron_end(src.intron_end),
        anchor_start(src.anchor_start), anchor_end(src.anchor_end),
        strand(src.strand), read_count(src.read_count) {
    }

    // Print BED with anchors as blocks and intron as gap.
    void print_anchor_bed(ostream& out) const {
        out << chrom
            << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << "jnc" << "\t" << read_count << "\t" << strand
            << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << "255,0,0" << "\t" << 2
            << "\t" << intron_start - anchor_start << "," << anchor_end - intron_end
            << "\t" << "0," << intron_end - anchor_start << endl;
    }

    // Print BED with intron as block
    void print_intron_bed(ostream& out) const {
        out << chrom
            << "\t" << intron_start << "\t" << intron_end
            << "\t" << "jnc" << "\t" << read_count << "\t" << strand << endl;
    }

};

// Compare two junctions
static inline bool compare_junctions(const Junction *j1,
                                     const Junction *j2) {
    if (j1->chrom < j2->chrom){
        return true;
    }
    if (j1->chrom > j2->chrom){
        return false;
    }
    // Same chromosome
    if (j1->anchor_start < j2->anchor_start) {
        return true;
    }
    if (j1->anchor_start > j2->anchor_start) {
        return false;
    }
    // Same start
    if (j1->anchor_end < j2->anchor_end) {
        return true;
    }
    if (j1->anchor_end > j2->anchor_end) {
        return false;
    }
    // Same end
    return j1->strand < j2->strand;
}

// The class that deals with creating the junctions
class JunctionsExtractor {
public:
    static const unsigned UNSTRANDED = 0;
    static const unsigned RF_STRANDED = 1;
    static const unsigned FR_STRANDED = 2;
private:
    // Minimum anchor length for junctions
    // Junctions need atleast this many bp overlap
    // on both ends.  Mismatch bases are not included in the count.
    uint32_t min_anchor_length_;
    // Minimum length of an intron, i.e min junction width
    uint32_t min_intron_length_;
    // Maximum length of an intron, i.e max junction width
    uint32_t max_intron_length_;
    //strandness of data; 0 = unstranded, 1 = RF, 2 = FR
    int strandness_;

    // Alignment file, used for error messsages
    string bam_;

    // target index to target (chrom) name
    vector<string> targets_;
    
    // Map to store the junctions by intron coordinates
    map<JunctionKey, Junction*> junctions_;

    // internal functions
    void save_targets(bam_hdr_t *header);
    bool junction_qc(uint32_t anchor_start, uint32_t intron_start,
                     uint32_t intron_end, uint32_t anchor_end,
                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt);
    int parse_alignment_into_junctions(bam1_t *aln);
    void add_junction(const string& chrom, char strand,
                      uint32_t anchor_start, uint32_t intron_start, 
                      uint32_t intron_end, uint32_t anchor_end);
    void process_junction(const string& chrom, char strand,
                          uint32_t anchor_start, uint32_t intron_start,
                          uint32_t intron_end, uint32_t anchor_end,
                          uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt);
    char get_junction_strand_XS(bam1_t *aln);
    char get_junction_strand_flag(bam1_t *aln);
    char get_junction_strand(bam1_t *aln);

    JunctionsExtractor() {
        assert(false); // Default constructor not allowed
    }
public:
    JunctionsExtractor(uint32_t min_anchor_length,
                       uint32_t min_intron_length,
                       uint32_t max_intron_length,
                       int strandness):
        min_anchor_length_(min_anchor_length),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        strandness_(strandness) {
    }
    
    // Identify exon-exon junctions
    void identify_junctions_from_bam(const string& bam);

    // Print BED with anchors as blocks and intron as gap.
    void print_anchor_bed(ostream& out);

    // Print BED with intron as block
    void print_intron_bed(ostream& out);

    // Get a vector of all the junctions
    vector<Junction*> get_junctions_sorted();

    // Get the BAM filename
    const string& get_bam() {
        return bam_;
    }
};

#endif
