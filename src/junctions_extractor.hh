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
#include <stdio.h>
#include <map>
#include <vector>
#include "htslib/sam.h"

using namespace std;

// Strandness of data; 0 = unstranded, 1 = RF, 2 = FR
// DO NOT CHANGE VALUES, code depends on it.
typedef enum {
    UNSTRANDED = 0,
    RF_STRANDED = 1,
    FR_STRANDED = 2
} Strandness;


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
    uint32_t anchor_start;  // This is the intron_start - max overhang
    uint32_t anchor_end;   // This is the end + max overhang
    char strand;
    // Number of reads supporting the junction, by category
    unsigned int read_counts[READ_CATEGORY_MAX + 1];

    Junction(const string& chrom1, uint32_t intron_start1, uint32_t intron_end1,
             uint32_t anchor_start1, uint32_t anchor_end1,
             char strand1):
        chrom(chrom1), intron_start(intron_start1), intron_end(intron_end1),
        anchor_start(anchor_start1), anchor_end(anchor_end1),
        strand(strand1)  {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = 0;
        }
    }
    Junction(const Junction& src):
        chrom(src.chrom), intron_start(src.intron_start), intron_end(src.intron_end),
        anchor_start(src.anchor_start), anchor_end(src.anchor_end),
        strand(src.strand)  {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = 0;
        }
    }

    // sum different read counts
    unsigned total_read_count() const {
        unsigned cnt = 0;
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            cnt += read_counts[i];
        }
        return cnt;
    }

    // Print BED with anchors as blocks and intron as gap.  ijunc is used to
    // make the BED name.
    void print_anchor_bed(unsigned ijunc,
                          ostream& out) const {

        if (chrom[0]>='0' && chrom[1]<='9')
                out << "chr";

        out << chrom
            << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << "sj" << ijunc << "\t" << total_read_count() << "\t" << strand
            << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << "255,0,0" << "\t" << 2
            << "\t" << intron_start - anchor_start << "," << anchor_end - intron_end
            << "\t" << "0," << intron_end - anchor_start << endl;
    }

    // Print BED with intron as block  ijunc is used to
    // make the BED name.
    void print_intron_bed(unsigned ijunc,
                          ostream& out) const {
        if (chrom[0]>='0' && chrom[1]<='9')
                out << "chr";

        out << chrom
            << "\t" << intron_start << "\t" << intron_end
            << "\t" << "sj" << ijunc << "\t" << total_read_count() << "\t" << strand << endl;
    }

    // Print header for junction call TSV
    static void print_juncion_call_header(ostream& out) {
        out << "chrom" << "\t" << "intron_start" << "\t" << "intron_end" << "\t" << "strand"
            << "\t" << "uniq_mapped_count" << "\t" << "multi_mapped_count" << "\t" << "unsure_mapped_count"
            << "\t" << "max_left_overhang" << "\t" << "max_right_overhang" << endl;
    }

    // Print row to junction call TSV
    void print_juncion_call_row(ostream& out) const {
        out << chrom << "\t" << intron_start << "\t" << intron_end << "\t" << strand
            << "\t" << read_counts[SINGLE_MAPPED_READ] << "\t" << read_counts[MULTI_MAPPED_READ] << "\t" << read_counts[UNSURE_READ]
            << "\t" << intron_start - anchor_start << "\t" << anchor_end - intron_end << endl;
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
    Strandness strandness_;

    // Alignment file
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
    void parse_alignment_into_junctions(bam1_t *aln);
    void process_alignment(bam1_t *aln);
    void add_junction(bam1_t *aln, const string& chrom, char strand,
                      uint32_t anchor_start, uint32_t intron_start, 
                      uint32_t intron_end, uint32_t anchor_end);
    void process_junction(bam1_t *aln, const string& chrom, char strand,
                          uint32_t anchor_start, uint32_t intron_start,
                          uint32_t intron_end, uint32_t anchor_end,
                          uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt);
    char get_junction_strand_XS(bam1_t *aln);
    char get_junction_strand_flag(bam1_t *aln);
    char get_junction_strand(bam1_t *aln);
    int get_num_aligns(bam1_t *aln);
    ReadCategory get_category_from_tag(bam1_t *aln);
    ReadCategory get_read_category(bam1_t *aln);
    samFile* open_pass_through(samFile *in_sam,
                               bam_hdr_t *in_header,
                               const string& bam_pass_through);

    JunctionsExtractor() {
        assert(false); // Default constructor not allowed
    }
public:
    JunctionsExtractor(uint32_t min_anchor_length,
                       uint32_t min_intron_length,
                       uint32_t max_intron_length,
                       Strandness strandness):
        min_anchor_length_(min_anchor_length),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        strandness_(strandness) {
    }
    
    // Identify exon-exon junctions
    void identify_junctions_from_bam(const string& bam,
                                     const string& bam_pass_through);

    // Print BED with anchors as blocks and intron as gap.
    void print_anchor_bed(ostream& out);

    // Print BED with intron as block
    void print_intron_bed(ostream& out);

    // Get a vector of all the junctions
    vector<Junction*> get_junctions_sorted();
};

#endif
