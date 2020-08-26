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

#ifndef JUNCTIONS_EXTRACTOR_H
#define JUNCTIONS_EXTRACTOR_H

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include "htslib/sam.h"
#include "junctions.hh"

class Genome;

using namespace std;

// Strandness of data; 0 = unstranded, 1 = RF, 2 = FR
// DO NOT CHANGE VALUES, code depends on it.
typedef enum {
    UNSTRANDED = 0,
    RF_STRANDED = 1,
    FR_STRANDED = 2
} Strandness;


// categories to exclude
typedef enum {
    EXCLUDE_NONE = 0,
    EXCLUDE_MULTI = 0x01,
} ExcludeCats;


// The class that deals with creating the junctions.
// Coordinate-sorted BAM are require to process one target
// at a time, reducing memory usage.
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

    // exclude categories
    unsigned excludes_;
    
    // Alignment file
    string bam_;

    // Object to get genome sequences, which maybe NULL
    Genome *genome_;

    // missing genomic targets
    bool skip_missing_targets_;
    
    // target index to target (chrom) name
    vector<string> targets_;
    
    // Map to store the junctions by intron coordinates
    JunctionTable junctions_;

    // missing genomic sequence that have been warned about
    set<string> missing_genomic_warned_;

    // update tags from observed orientation
    bool set_XS_strand_tag_;
    bool set_TS_strand_tag_;

    // checking for not sorted
    set<int> done_targets_;

    // alignment buffer, used to handle one level of read-ahead
    bam1_t *aln_buf_;
    bool aln_pending_;

    // used for generate BED names
    uint32_t ijunc_;

    // BAM files
    samFile *in_sam_;
    bam_hdr_t *in_header_;
    samFile *out_sam_;
    
    // debugging trace output if not NULL
    ostream *trace_fh_;
    
    // internal functions
    void save_targets(bam_hdr_t *header);
    bool junction_qc(bam1_t *aln, uint32_t anchor_start, uint32_t intron_start,
                     uint32_t intron_end, uint32_t anchor_end,
                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt);
    void parse_alignment_into_junctions(bam1_t *aln,
                                        int *orientCnt);
    void process_alignment(bam1_t *aln);
    void write_pass_through(bam1_t *aln,
                            bam_hdr_t *in_header,
                            samFile* out_sam);
    
    Junction *create_junction(bam1_t *aln, const JunctionKey &key,
                              const string& chrom, char strand,
                              uint32_t anchor_start, uint32_t intron_start,
                              uint32_t intron_end, uint32_t anchor_end);
    Junction *update_junction(bam1_t *aln, const JunctionKey &key,
                              const string& chrom, char strand,
                              uint32_t anchor_start, uint32_t intron_start,
                              uint32_t intron_end, uint32_t anchor_end);
    Junction *add_junction(bam1_t *aln, const string& chrom, char strand,
                           uint32_t anchor_start, uint32_t intron_start, 
                           uint32_t intron_end, uint32_t anchor_end);
    void update_strand_tag(const char* tag, char valType, int orientCnt, bam1_t *aln);
    void process_junction(bam1_t *aln, const string& chrom, char strand,
                          uint32_t anchor_start, uint32_t intron_start,
                          uint32_t intron_end, uint32_t anchor_end,
                          uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt,
                          int *orientCnt);
    char get_junction_strand_XS(bam1_t *aln);
    char get_junction_strand_flag(bam1_t *aln);
    char get_junction_strand(bam1_t *aln);
    int get_num_aligns(bam1_t *aln);
    ReadCategory get_category_from_tag(bam1_t *aln);
    ReadCategory get_read_category(bam1_t *aln);
    samFile* open_pass_through(samFile *in_sam,
                               bam_hdr_t *in_header,
                               const string& bam_pass_through);
    bool is_canonical(const string& junctions);
    string get_splice_sites(const string& chrom, char strand,
                            uint32_t intron_start, uint32_t intron_end);

    bam1_t* read_align();
    void check_new_target(bam1_t *aln);

    JunctionsExtractor() {
        assert(false); // Default constructor not allowed
    }
public:
    JunctionsExtractor(uint32_t min_anchor_length,
                       uint32_t min_intron_length,
                       uint32_t max_intron_length,
                       Strandness strandness,
                       unsigned excludes,
                       Genome *genome,
                       bool skip_missing_targets,
                       bool set_XS_strand_tag,
                       bool set_TS_strand_tag,
                       ostream *trace_fh):
        min_anchor_length_(min_anchor_length),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        strandness_(strandness),
        excludes_(excludes),
        genome_(genome),
        skip_missing_targets_(skip_missing_targets),
        set_XS_strand_tag_(set_XS_strand_tag),
        set_TS_strand_tag_(set_TS_strand_tag),
        aln_buf_(NULL),
        aln_pending_(false),
        ijunc_(0),
        in_sam_(NULL),
        in_header_(NULL),
        out_sam_(NULL),
        trace_fh_(trace_fh) {
    }

    ~JunctionsExtractor();

    // free junctions found so far
    void clear() {
        junctions_.clear();
    }
    
    // Open BAMs
    void open(const string& bam,
              const string& bam_pass_through);

    // get the number of targets
    int get_num_targets() const {
        return targets_.size();
    }
    
    // Identify exon-exon junctions
    void identify_junctions_for_target(int target_index);

    // Copy remaining reads (unaligned) to pass-through, if it is open
    void copy_unaligned_reads();

    // Get a vector of all the junctions
    JunctionVector get_junctions() {
        return junctions_.get_junctions();
    }
};

// open a output file if not empty, otherwise return NULL
inline ostream* open_out_or_null(const string& fname) {
    return (fname != "") ? new ofstream(fname.c_str()) : NULL;
}

#endif
