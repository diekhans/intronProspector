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

// Data save for an intron
class Junction {
public:
    const string& chrom;
    uint32_t start;
    uint32_t end;
    char strand;
    // Number of reads supporting the junction
    unsigned int read_count;
    // This is the start - max overhang
    uint32_t thick_start;
    // This is the end + max overhang
    uint32_t thick_end;

    Junction(const string& chrom1, uint32_t start1, uint32_t end1,
             uint32_t thick_start1, uint32_t thick_end1,
             char strand1):
        chrom(chrom1), start(start1), end(end1),
        thick_start(thick_start1), thick_end(thick_end1),
        strand(strand1), read_count(0) {
    }
    Junction(const Junction& src):
        chrom(src.chrom), start(src.start), end(src.end),
        thick_start(src.thick_start), thick_end(src.thick_end),
        strand(src.strand), read_count(src.read_count) {
    }

    // Print junction
    void print(ostream& out) const {
        out << chrom <<
            "\t" << thick_start << "\t" << thick_end <<
            "\t" << "jnc" << "\t" << read_count << "\t" << strand <<
            "\t" << thick_start << "\t" << thick_end <<
            "\t" << "255,0,0" << "\t" << 2 <<
            "\t" << start - thick_start << "," << thick_end - end <<
            "\t" << "0," << end - thick_start << endl;
    }

};

// Compare two junctions
// Return true if j1.start < j2.start
// If j1.start == j2.start, return true if j1.end < j2.end
static inline bool compare_junctions(const Junction &j1,
                                     const Junction &j2) {
    // Different chromosome
    if (j1.chrom < j2.chrom){
        return true;
    }
    if (j1.chrom > j2.chrom){
        return false;
    }
    // Same chromosome
    if (j1.thick_start < j2.thick_start) {
        return true;
    }
    if (j1.thick_start > j2.thick_start) {
        return false;
    }
    if (j1.thick_end < j2.thick_end) {
        return true;
    }
    return (j1.thick_end > j2.thick_end);
}

// Sort a vector of junctions
template <class CollectionType>
inline void sort_junctions(CollectionType &junctions) {
    sort(junctions.begin(), junctions.end(), compare_junctions);
}

// The class that deals with creating the junctions
class JunctionsExtractor {
public:
    static const uint32_t DEFAULT_MIN_ANCHOR_LENGTH = 8;
    static const uint32_t DEFAULT_MIN_INTRON_LENGTH = 70;
    static const uint32_t DEFAULT_MAX_INTRON_LENGTH = 500000;

    static const int UNSTRANDED = 0;
    static const int RF_STRANDED = 1;
    static const int FR_STRANDED = 2;
private:
    // Alignment file
    string bam_;
    // Minimum anchor length for junctions
    // Junctions need atleast this many bp overlap
    // on both ends.
    uint32_t min_anchor_length_;
    // Minimum length of an intron, i.e min junction width
    uint32_t min_intron_length_;
    // Maximum length of an intron, i.e max junction width
    uint32_t max_intron_length_;
    //strandness of data; 0 = unstranded, 1 = RF, 2 = FR
    int strandness_;

    // target index to target (chrom) name
    vector<string> targets_;
    
    // Map to store the junctions
    // The key is "chr:start-end:strand"
    // The value is an object of type Junction(see above)
    map<string, Junction*> junctions_;
    // Maintain a sorted list of junctions
    //FIXME: vector<Junction> junctions_vector_;

    // internal functions
    void save_targets(bam_hdr_t *header);
    bool junction_qc(uint32_t anchor_start, uint32_t anchor_end,
                     uint32_t thick_start, uint32_t thick_end);
    string make_junction_key(const string& chrom, char strand,
                             uint32_t start, uint32_t end);
    int parse_alignment_into_junctions(bam_hdr_t *header, bam1_t *aln);
    int parse_read(bam_hdr_t *header, bam1_t *aln);
    int parse_cigar_into_junctions(string chr, int read_pos,
                                   uint32_t *cigar, int n_cigar);
    void add_junction(const string& chrom, char strand,
                      uint32_t anchor_start, uint32_t anchor_end,
                      uint32_t thick_start, uint32_t thick_end);
    char get_junction_strand_XS(bam1_t *aln);
    char get_junction_strand_flag(bam1_t *aln);
    char get_junction_strand(bam1_t *aln);
    void create_junctions_vector();
    
    JunctionsExtractor() {
        assert(false); // Default constructor not allowed
    }
public:
    JunctionsExtractor(string bam,
                       uint32_t min_anchor_length=DEFAULT_MIN_ANCHOR_LENGTH,
                       uint32_t min_intron_length=DEFAULT_MIN_INTRON_LENGTH,
                       uint32_t max_intron_length=DEFAULT_MAX_INTRON_LENGTH,
                       int strandness = RF_STRANDED):
        bam_(bam),
        min_anchor_length_(min_anchor_length),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        strandness_(strandness) {
    }
    
    // Identify exon-exon junctions
    void identify_junctions_from_bam();

    // Print all the junctions
    void print_all_junctions(ostream& out=cout);

    // Get a vector of all the junctions
    vector<Junction> get_all_junctions();

    // Get the BAM filename
    const string& get_bam() {
        return bam_;
    }
};

#endif
