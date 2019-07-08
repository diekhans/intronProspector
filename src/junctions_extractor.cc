/*  junctions_extractor.cc -- splice junction extraction.

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

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "junctions_extractor.hh"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

using namespace std;

float Junction::NULL_CONFIDENCE = nanf("not-a-number");

// convert a char to a string
static string char_to_string(char ch) {
    char chs [2] = {ch, '\0'};
    return string(chs);
}

#if 1
// convert an integer to a string
static string int_to_string(int num) {
    stringstream s1;
    s1 << num;
    return s1.str();
}
#endif

// lazy calculation of confidence
void Junction::lazy_get_confidence() {
    assert(isnan(confidence));
    if (total_read_count() == 0) {
        confidence = 0.0;
    } else {
        confidence = calculate_confidence();
    }
}

// calculation of confidence
float Junction::calculate_confidence() {
    // put in order for counting
    std::sort(read_offsets.begin(), read_offsets.end());
    float summation = 0.0;
    for (unsigned i = 0; i < read_offsets.size(); ) {
        uint16_t pos = read_offsets[i];
        unsigned cnt = 0;
        for (; (i < read_offsets.size()) && (read_offsets[i] == pos); i++) {
            cnt++;
        }
        float p = float(cnt) / read_offsets.size();
        if (p != 0.0) {
            summation += p * log2f(p);
        }
    }
    return (summation == 0.0) ? 0.0 : -summation;
}

// Sort a vector of junctions
template <class CollectionType>
static inline void sort_junctions(CollectionType &junctions) {
    sort(junctions.begin(), junctions.end(), compare_junctions);
}

// Do some basic qc on the junction
bool JunctionsExtractor::junction_qc(uint32_t anchor_start, uint32_t intron_start,
                                     uint32_t intron_end, uint32_t anchor_end,
                                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt) {
    if ((intron_end - intron_start < min_intron_length_) ||
        (intron_end - intron_start > max_intron_length_)) {
        return false;
    } else if (((intron_start - anchor_start) - left_mismatch_cnt < min_anchor_length_) ||
               ((anchor_end - intron_end) - right_mismatch_cnt < min_anchor_length_)) {
        return false;
    } else {
        return true;
    }
}

// Add a junction to the junctions map
void JunctionsExtractor::add_junction(bam1_t *aln, const string& chrom, char strand,
                                      uint32_t anchor_start, uint32_t intron_start,
                                      uint32_t intron_end, uint32_t anchor_end) {
    // Construct key chr:start-end:strand
    JunctionKey key(chrom, intron_start, intron_end, strand);

    // Check if new junction
    Junction *junc = NULL;
    if (!junctions_.count(key)) {
        junc = new Junction(chrom, intron_start, intron_end, anchor_start, anchor_end, strand);
        junctions_[key] = junc;
    } else {
         // existing junction
        junc = junctions_[key];
        // Check if thick starts are any better
        if (anchor_start < junc->anchor_start)
            junc->anchor_start = anchor_start;
        if (anchor_end > junc->anchor_end)
            junc->anchor_end = anchor_end;
    }
    junc->count_read(get_read_category(aln),
                     intron_start - aln->core.pos);
}

// Validate a junction and save if it passes.
void JunctionsExtractor::process_junction(bam1_t *aln, const string& chrom, char strand,
                                          uint32_t anchor_start, uint32_t intron_start,
                                          uint32_t intron_end, uint32_t anchor_end,
                                          uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt) {
    if (junction_qc(anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt)) {
        add_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end);
    }
}

// Print all the junctions - this function needs work
vector<Junction*> JunctionsExtractor::get_junctions_sorted() {
    vector<Junction*> juncs;
    for (map<JunctionKey, Junction*>::iterator it = junctions_.begin(); it != junctions_.end(); it++) {
        juncs.push_back(it->second);
    }
    sort_junctions(juncs);
    return juncs;
}

// Get the strand from the XS aux tag
char JunctionsExtractor::get_junction_strand_XS(bam1_t *aln) {
    uint8_t *p = bam_aux_get(aln, "XS");
    if (p != NULL) {
        char strand = bam_aux2A(p);
        return strand ? strand : '.';
    } else {
        return '.';
    }
}

// Get the strand from the bitwise flag
char JunctionsExtractor::get_junction_strand_flag(bam1_t *aln) {
    uint32_t flag = (aln->core).flag;
    int reversed = bam_is_rev(aln);
    int mate_reversed = bam_is_mrev(aln);
    int first_in_pair = (flag & BAM_FREAD1) != 0;
    int second_in_pair = (flag & BAM_FREAD2) != 0;
    // strandness_ is 0 for unstranded, 1 for RF, and 2 for FR
    int bool_strandness = strandness_ - 1;
    int first_strand = !bool_strandness ^ first_in_pair ^ reversed;
    int second_strand = !bool_strandness ^ second_in_pair ^ mate_reversed;
    char strand;
    if (first_strand) {
        strand = '+';
    } else {
        strand = '-';
    }
    // if strand inferences from first and second in pair don't agree, we've got a problem
    if (first_strand == second_strand) {
        return strand;
    } else {
        return '.';
    }
}

// Get the strand
char JunctionsExtractor::get_junction_strand(bam1_t *aln) {
    if (strandness_ != UNSTRANDED) {
        return get_junction_strand_flag(aln);
    } else {
        return get_junction_strand_XS(aln);
    }
}

// get number of alignments for current read (NH tag), or -1 if not set,
// and 0 if there is some other issue (don't check error, just ignored).
int JunctionsExtractor::get_num_aligns(bam1_t *aln) {
    uint8_t *tag = bam_aux_get(aln, "NH");
    if (tag == NULL) {
        return -1;
    } else {
        return bam_aux2i(tag);
    }
}

// Determine category base on tag
ReadCategory JunctionsExtractor::get_category_from_tag(bam1_t *aln) {
    int num_aligns = get_num_aligns(aln);
    if (num_aligns <= 0) {
        return UNSURE_READ;
    } else if (num_aligns == 1) {
        return SINGLE_MAPPED_READ;
    } else {
        return MULTI_MAPPED_READ;
    }
}


// Determine the category of read.  Logic from samtools flagstat
ReadCategory JunctionsExtractor::get_read_category(bam1_t *aln) {
    unsigned flag = aln->core.flag;
    if (flag & (BAM_FQCFAIL | BAM_FUNMAP | BAM_FDUP)) {
        throw logic_error("should never look at FQCAIL, FUNMAP, FDUP reads");
    } else if (flag & BAM_FSECONDARY) {
        return MULTI_MAPPED_READ;  // secondary alignment
    } else if (flag & BAM_FSUPPLEMENTARY ) {
        return MULTI_MAPPED_READ;  // supplementary alignment
    } else {
        if (flag & BAM_FPAIRED) {
            // paired reads
            if ((flag & BAM_FPROPER_PAIR) && !(flag & BAM_FUNMAP)) {
                return get_category_from_tag(aln); // proper pair
            } else if (flag & BAM_FMUNMAP) {
                return UNSURE_READ; // mate unmapped
            } else if (aln->core.mtid != aln->core.tid) {
                return UNSURE_READ; // mapped to different chromosomes
            } else {
                return UNSURE_READ; // mapped to same chromosome, but not a proper pair
            }
        } else {
            return get_category_from_tag(aln);  // mapped, unpaired read
        }
    }
}

// Parse junctions from the read and store in junction map
void JunctionsExtractor::parse_alignment_into_junctions(bam1_t *aln) {
    const string& chrom = targets_[aln->core.tid];
    char strand = get_junction_strand(aln);
    int read_pos = aln->core.pos;
    int n_cigar = aln->core.n_cigar;
    uint32_t *cigar = bam_get_cigar(aln);
    uint32_t anchor_start = read_pos;
    uint32_t intron_start = read_pos;
    uint32_t intron_end = 0;
    uint32_t anchor_end = 0;
    uint32_t left_mismatch_cnt = 0;   // mismatches in left anchor
    uint32_t right_mismatch_cnt = 0;  // mismatches in right anchor
    bool started_junction = false;
    for (int i = 0; i < n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        switch(op) {
            case 'N': // skipped region from the reference
                if (!started_junction) {
                    intron_end = intron_start + len;
                    anchor_end = intron_end;
                    // Start the first one and remains started
                    started_junction = true;
                } else {
                    // Add the previous junction
                    process_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt);
                    anchor_start = intron_end;
                    intron_start = anchor_end;
                    intron_end = intron_start + len;
                    anchor_end = intron_end;
                    // For clarity - the next junction is now open
                    started_junction = true;
                }
                break;
            case '=':  // sequence match
            case 'M':  // alignment match (can be a sequence match or mismatch)
                if (!started_junction) {
                    intron_start += len;
                } else {
                    anchor_end += len;
                }
                break;
            case 'X':  // sequence mismatch
                if (!started_junction) {
                    intron_start += len;
                    left_mismatch_cnt += len;
                } else {
                    right_mismatch_cnt += len;
                    anchor_end += len;
                }
                break;
            case 'D':  // deletion from the reference
                if (!started_junction) {
                    intron_start += len;
                    anchor_start = intron_start;
                } else {
                    process_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt);
                    // Don't include these in the next anchor
                    intron_start = anchor_end + len;
                    anchor_start = intron_start;
                }
                started_junction = false;
                break;
            case 'I':  // insertion to the reference
            case 'S':  // soft clipping (clipped sequences present in SEQ)
                if (!started_junction) {
                    anchor_start = intron_start;
                } else {
                    process_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt);
                    // Don't include these in the next anchor
                    intron_start = anchor_end;
                    anchor_start = intron_start;
                }
                started_junction = false;
                break;
            case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                break;
            default:
                throw std::invalid_argument("Unknown cigar operation '" + char_to_string(op) + "' found in " + bam_);
        }
    }
    if (started_junction) {
        process_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt);
    }
}

// Process an alignment
void JunctionsExtractor::process_alignment(bam1_t *aln) {
    // skip if unmapped, only one cigar operation exists (likely all matches),
    // low-quality or duplicate
    if (!((aln->core.n_cigar <= 1) || (aln->core.tid < 0) || (aln->core.flag & (BAM_FQCFAIL | BAM_FDUP)))) {
        parse_alignment_into_junctions(aln);
    }
}


// build target array from bam header
void JunctionsExtractor::save_targets(bam_hdr_t *header) {
    for (int i = 0; i < header->n_targets; i++) {
        targets_.insert(targets_.end(), string(header->target_name[i]));
    }
}

// open pass-through file
samFile* JunctionsExtractor::open_pass_through(samFile *in_sam,
                                               bam_hdr_t *in_header,
                                               const string& bam_pass_through) {
    const htsFormat *fmt = hts_get_format(in_sam);
    if (!((fmt->format == sam) || (fmt->format == bam))) {
        throw invalid_argument("Error: pass-through is only implemented for SAM or BAM files: " + bam_pass_through);
    }
        
    samFile *out_sam = hts_open_format(bam_pass_through.c_str(), "w", fmt);
    if (out_sam == NULL) {
        throw runtime_error("Error opening BAM/SAM/CRAM file: " + bam_pass_through);
    }
    if (sam_hdr_write(out_sam, in_header) < 0) {
        throw runtime_error("Error writing SAM header: " + bam_pass_through);
    }
    return out_sam;
}

// The workhorse - identifies junctions from BAM
void JunctionsExtractor::identify_junctions_from_bam(const string& bam,
                                                     const string& bam_pass_through) {
    bam_ = bam;
    samFile *in_sam = sam_open(bam_.c_str(), "r");
    if (in_sam == NULL) {
        throw runtime_error("Error opening BAM/SAM/CRAM file: " + bam_);
    }
    bam_hdr_t *in_header = sam_hdr_read(in_sam);
    save_targets(in_header);

    samFile* out_sam = NULL;
    if (bam_pass_through != "") {
        out_sam = open_pass_through(in_sam, in_header, bam_pass_through);
    }

    bam1_t *aln = bam_init1();
    int stat;
    while((stat = sam_read1(in_sam, in_header, aln)) >= 0) {
        process_alignment(aln);
        if (out_sam != NULL) {
            if (sam_write1(out_sam, in_header, aln) < 0) {
                throw runtime_error("Error writing BAM record: " + bam_pass_through);
            }
        }
    }
    if (stat < -1) {
        throw runtime_error("Error reader BAM record: " + bam);
    }
    bam_destroy1(aln);
    if (out_sam != NULL) {
        sam_close(out_sam);
    }
    bam_hdr_destroy(in_header);
    sam_close(in_sam);
}
