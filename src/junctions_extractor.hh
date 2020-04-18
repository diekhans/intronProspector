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
#include <sstream>
#include <stdio.h>
#include <map>
#include <vector>
#include <math.h>
#include "htslib/sam.h"
class Genome;

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

#if 0 // UNUSED
// index of strand in arrays
typedef enum {
    UNKNOWN_IDX = 0,
    PLUS_IDX = 1,
    MINUS_IDX = 2,
    NUM_STRANDS = 3,
} StrandIdx;

static inline StrandIdx strand_to_idx(char strand) {
    switch (strand) {
        case '+':
            return PLUS_IDX;
        case '-':
            return MINUS_IDX;
        case '.':
            return UNKNOWN_IDX;
        default:
            assert(false);
            return UNKNOWN_IDX;
    }
}

static inline char idx_to_strand(StrandIdx strand_idx) {
    switch (strand_idx) {
        case PLUS_IDX:
            return '+';
        case MINUS_IDX:
            return '+';
        case UNKNOWN_IDX:
            return '.';
        default:
            assert(false);
            return '.';
    }
}
#endif

// categories to exclude
typedef enum {
    EXCLUDE_NONE = 0,
    EXCLUDE_MULTI = 0x01,
} ExcludeCats;


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
    
    const string& chrom;
    uint32_t intron_start;
    uint32_t intron_end;
    uint32_t anchor_start;  // This is the intron_start - max overhang
    uint32_t anchor_end;    // This is the end + max overhang
    char strand;
    string splice_sites;    // splite sites in the form GT/AG if genome available

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
             char strand1, const string& splice_sites):
        chrom(chrom1), intron_start(intron_start1), intron_end(intron_end1),
        anchor_start(anchor_start1), anchor_end(anchor_end1),
        strand(strand1), confidence(NULL_CONFIDENCE), splice_sites(splice_sites) {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = 0;
        }
    }
    Junction(const Junction& src):
        chrom(src.chrom), intron_start(src.intron_start), intron_end(src.intron_end),
        anchor_start(src.anchor_start), anchor_end(src.anchor_end),
        strand(src.strand), read_offsets(src.read_offsets), confidence(src.confidence),
        splice_sites(splice_sites)  {
        for (unsigned i = 0; i <= READ_CATEGORY_MAX; i++) {
            read_counts[i] = 0;
        }
    }

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
    
    // Trim read counts to fit in BED score restriction or 0..1000
    static unsigned read_count_to_bed_score(uint64_t read_count) {
        return (read_count <= 1000) ? read_count : 1000;
    }

    // naive mapping to ucsc chrom name
    const string make_ucsc_chrom(const string& chrom) const {
        if ((chrom.find_first_not_of("0123456789") == std::string::npos)
            || (chrom == "X") || (chrom == "Y")) {
            return string("chr") + chrom;
        } else if (chrom == "MT") {
            return "chrM";
        } else {
            return chrom;
        }
    }

    // color-code to use for splice-sites
    const string& getBedColor() const {
        static const string NO_INFO_COLOR = "64,128,64";
        static const string U2_COLOR = "0,128,0";
        static const string U12_COLOR = "0,0,128";
        static const string UNKNOWN_COLOR = "128,0,0";
        if (splice_sites == "") {
            return NO_INFO_COLOR;
        }
        if (splice_sites == "GT/AG") {
            return U2_COLOR;
        }
        if ((splice_sites == "AT/AC") or
            (splice_sites == "AT/AG") or 
            (splice_sites == "GT/AC") or
            (splice_sites == "GT/AG")) {
            return U12_COLOR;
        }
        return UNKNOWN_COLOR;
    }
    
    // Print BED with anchors as blocks and intron as gap.  ijunc is used to
    // make the BED name.
    void print_anchor_bed(unsigned ijunc,
                          bool map_to_ucsc,
                          ostream& out) const {
        if (map_to_ucsc) {
            out << make_ucsc_chrom(chrom);
        } else {
            out << chrom;
        }        
        out << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << "sj" << ijunc;
        if (splice_sites != "") {
            out << '_' << splice_sites;
        }
        out << "\t" << read_count_to_bed_score(total_read_count()) << "\t" << strand
            << "\t" << anchor_start << "\t" << anchor_end
            << "\t" << getBedColor() << "\t" << 2
            << "\t" << intron_start - anchor_start << "," << anchor_end - intron_end
            << "\t" << "0," << intron_end - anchor_start << endl;
    }

    // Print BED with intron as block  ijunc is used to
    // make the BED name.
    void print_intron_bed(unsigned ijunc,
                          bool map_to_ucsc,
                          ostream& out) const {
        if (map_to_ucsc) {
            out << make_ucsc_chrom(chrom);
        } else {
            out << chrom;
        }        
        out << "\t" << intron_start << "\t" << intron_end
            << "\t" << "sj" << ijunc;
        if (splice_sites != "") {
            out << '_' << splice_sites;
        }
        out << "\t" << read_count_to_bed_score(total_read_count()) << "\t" << strand
            << "\t" << intron_start << "\t" << intron_end << "\t" << getBedColor() << endl;
    }

    // Print header for junction call TSV
    static void print_junction_call_header(bool have_genome,ostream& out) {
        out << "chrom" << "\t" << "intron_start" << "\t" << "intron_end" << "\t" << "strand"
            << "\t" << "uniq_mapped_count" << "\t" << "multi_mapped_count" << "\t" << "unsure_mapped_count"
            << "\t" << "max_left_overhang" << "\t" << "max_right_overhang" << "\t" << "confidence";
        if (have_genome) {
            out << "\t" << "splice_sites";
        }
        out << endl;
    }

    // Print row to junction call TSV
    void print_junction_call_row(bool have_genome,
                                 bool map_to_ucsc,
                                 ostream& out) const {
        if (map_to_ucsc) {
            out << make_ucsc_chrom(chrom);
        } else {
            out << chrom;
        }        
        out << "\t" << intron_start << "\t" << intron_end << "\t" << strand
            << "\t" << read_counts[SINGLE_MAPPED_READ] << "\t" << read_counts[MULTI_MAPPED_READ] << "\t" << read_counts[UNSURE_READ]
            << "\t" << intron_start - anchor_start << "\t" << anchor_end - intron_end
            << "\t" << get_confidence();
        if (have_genome) {
            out << "\t" << splice_sites;
        }
        out << endl;
    }

};

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
    map<JunctionKey, Junction*> junctions_;

    // missing genomic sequence that have been warned about
    map<string, bool> missing_genomic_warned_;

    // debugging trace output if not NULL
    ostream *trace_fh_;
    
    // internal functions
    void save_targets(bam_hdr_t *header);
    bool junction_qc(bam1_t *aln, uint32_t anchor_start, uint32_t intron_start,
                     uint32_t intron_end, uint32_t anchor_end,
                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt);
    void parse_alignment_into_junctions(bam1_t *aln);
    void process_alignment(bam1_t *aln);
    void create_junction(bam1_t *aln, const JunctionKey &key,
                         const string& chrom, char strand,
                         uint32_t anchor_start, uint32_t intron_start,
                         uint32_t intron_end, uint32_t anchor_end);
    void update_junction(bam1_t *aln, const JunctionKey &key,
                         const string& chrom, char strand,
                         uint32_t anchor_start, uint32_t intron_start,
                         uint32_t intron_end, uint32_t anchor_end);
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
    bool is_canonical(const string& junctions);
    string get_splice_sites(const string& chrom, char strand,
                            uint32_t intron_start, uint32_t intron_end);

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
                       ostream *trace_fh):
        genome_(genome),
        skip_missing_targets_(skip_missing_targets),
        min_anchor_length_(min_anchor_length),
        min_intron_length_(min_intron_length),
        max_intron_length_(max_intron_length),
        strandness_(strandness),
        excludes_(excludes),
        trace_fh_(trace_fh) {
    }
    
    // Identify exon-exon junctions
    void identify_junctions_from_bam(const string& bam,
                                     const string& bam_pass_through);

    // Print BED with anchors as blocks and intron as gap.
    void print_anchor_bed(ostream& out);

    // Print BED with intron as block
    void print_intron_bed(ostream& out);

    // Get a vector of all the junctions
    vector<Junction*> get_junctions();
};

// sort by anchor start
void sort_by_anchors(vector<Junction*>& junctions);

// sort by intron start
void sort_by_introns(vector<Junction*>& junctions);

#endif
