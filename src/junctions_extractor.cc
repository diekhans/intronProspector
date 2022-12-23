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
#include <stdexcept>
#include <set>
#include <algorithm>
#include "junctions_extractor.hh"
#include "type_ops.hh"
#include "genome.hh"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

using namespace std;

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

static bool is_canonical(const string& junctions) {
    static set<string> canonicals;  // has all combinations of upper and lower case
    if (canonicals.size() == 0) {
        build_canonical(canonicals);
    }
    return canonicals.count(junctions) > 0;
}

// get number of alignments for current read (NH tag), or -1 if not set,
// and 0 if there is some other issue (don't check error, just ignored).
static int get_num_aligns(bam1_t *aln) {
    uint8_t *tag = bam_aux_get(aln, "NH");
    if (tag == NULL) {
        return -1;
    } else {
        return bam_aux2i(tag);
    }
}

// Determine category base on tag
static ReadCategory get_category_from_tag(bam1_t *aln) {
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
static ReadCategory get_read_category(bam1_t *aln) {
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

// update the specific strand tag from the observed orientation
static void update_strand_tag(const char* tag, char valType, int orientCnt, bam1_t *aln) {
    const uint8_t strand = (orientCnt > 0) ? '+' : '-';
    uint8_t *existing = bam_aux_get(aln, tag);
    if (existing != NULL) {
        int err = bam_aux_del(aln, existing);
        if (err != 0) {
            throw runtime_error("Error deleting read tag");
        }
    }
    int err = bam_aux_append(aln, tag, valType, 1, &strand);
    if (err != 0) {
        throw runtime_error("Error adding read tag");
    }
}


typedef enum {
    IN_LIMBO = 0x1,
    IN_LEFT_ANCHOR = 0x2,
    IN_INTRON = 0x4,
    IN_RIGHT_ANCHOR = 0x8
} State;

/*
 * state machine for calling introns from a read.
 */
class ReadJunctionExtractor {
    private:
    bam1_t *aln_;
    const string& chrom_;
    char strand_;
    JunctionsExtractor *junc_extractor_;
    State state_;
    bool started_junction_;
    hts_pos_t target_pos_;
    hts_pos_t anchor_start_;   // -1 when not in anchor/intron
    hts_pos_t intron_start_;
    hts_pos_t intron_end_;
    hts_pos_t anchor_end_;
    uint32_t left_mismatch_cnt_;   // mismatches in anchors
    uint32_t right_mismatch_cnt_;
    uint32_t left_indel_cnt_;   // indels in anchors
    uint32_t right_indel_cnt_;
    uint32_t left_indel_max_;   // max length of indels in anchors
    uint32_t right_indel_max_;
    int orientCnt_;  // sum of +1 for positive strand introns, -1 for negative strand

    void record_junction() {
        junc_extractor_->record_junction(aln_, chrom_, strand_,
                                         anchor_start_, intron_start_, intron_end_, anchor_end_,
                                         left_mismatch_cnt_, right_mismatch_cnt_,
                                         left_indel_cnt_, right_indel_cnt_,
                                         left_indel_max_, right_indel_max_, &orientCnt_);
    }

    void enter_limbo() {
        if (state_ == IN_RIGHT_ANCHOR) {
            record_junction();
        }
        state_ = IN_LIMBO;
    }

    void shift_anchor_counts() {
        left_mismatch_cnt_ = right_mismatch_cnt_;
        right_mismatch_cnt_ = 0;
        left_indel_cnt_ = right_indel_cnt_;
        right_indel_cnt_ = 0;
        left_indel_max_ = right_indel_max_;
        right_indel_max_ = 0;
    }
    
    void count_mismatchs(uint32_t len) {
        if (state_ == IN_LEFT_ANCHOR) {
            left_mismatch_cnt_ += len;
        } else if (state_ == IN_RIGHT_ANCHOR) {
            right_mismatch_cnt_ += len;
        }
    }
    
    void count_indels(uint32_t len) {
        if (state_ == IN_LEFT_ANCHOR) {
            left_indel_cnt_ += len;
            left_indel_max_ = max(left_indel_max_, len);
        } else if (state_ == IN_RIGHT_ANCHOR) {
            right_indel_cnt_ += len;
            right_indel_max_ = max(right_indel_max_, len);
        }
    }
    
    public:
    ReadJunctionExtractor(bam1_t *aln,
                          const string& chrom,
                          JunctionsExtractor *junc_extractor):
        aln_(aln),
        chrom_(chrom),
        strand_(junc_extractor->get_junction_strand(aln_)),
        junc_extractor_(junc_extractor),
        state_(IN_LIMBO),
        target_pos_(aln->core.pos),
        anchor_start_(-1),
        intron_start_(-1),
        intron_end_(-1),
        anchor_end_(-1),
        left_mismatch_cnt_(0),
        right_mismatch_cnt_(0),
        left_indel_cnt_(0),
        right_indel_cnt_(0),
        left_indel_max_(0),
        right_indel_max_(0),
        orientCnt_(0) {
    }

    // op: N
    void enter_intron(char op, uint32_t len) {
        //if in right anchor, close previous intron, move right anchor to left anchor
        //if in left anchor, record intron_range
        //elfi ignore

        if (state_ == IN_RIGHT_ANCHOR) {
            // right anchor of previous exon is left anchor of new one,
            record_junction();
            anchor_start_ = intron_end_;  // previous intron
            shift_anchor_counts();
            state_ = IN_LEFT_ANCHOR; 
        }

        // this will ignore intron if there is not a left anchor
        if (state_ == IN_LEFT_ANCHOR) {
            intron_start_ = target_pos_;
            intron_end_ = target_pos_ + len;
            state_ = IN_INTRON;
        }
        target_pos_ += len;
    }

    // op: =,M,X
    void enter_aligned(char op, uint32_t len) {
        // aligned region, becomes part of an anchor
        if (state_ == IN_LIMBO) {
            // start new left-anchor
            anchor_start_ = target_pos_;
            anchor_end_ = target_pos_ + len;
            left_mismatch_cnt_ = right_mismatch_cnt_ = 0;
            state_ = IN_LEFT_ANCHOR;
        } else if (state_ & (IN_LEFT_ANCHOR | IN_RIGHT_ANCHOR)) {
            // extend left or right anchors
            anchor_end_ = target_pos_ + len;
        } else if (state_ == IN_INTRON) {
            // intron ended, starting right anchor
            intron_end_ = target_pos_;
            anchor_end_= target_pos_ + len;
            state_ = IN_RIGHT_ANCHOR;
        }
        if (op == 'X') {
            count_mismatchs(len);
        }
        target_pos_ += len;
    }
    
    // op: D
    void enter_query_delete(char op, uint32_t len) {
        // non-intron insertion in genome
        count_indels(len);
        enter_limbo();
        target_pos_ += len;
    }
    
    // op: I
    void enter_query_insert(char op, uint32_t len) {
        // query insertion
        count_indels(len);
        enter_limbo();
    }

    // end of read
    void finish() {
        if (state_ == IN_RIGHT_ANCHOR) {
            record_junction();
        }
    }
    
    /* parse cigar string into junctions */
    void parse_read(int *orientCnt) {
        int n_cigar = aln_->core.n_cigar;
        uint32_t *cigar = bam_get_cigar(aln_);
        for (int i = 0; i < n_cigar; ++i) {
            char op = bam_cigar_opchr(cigar[i]);
            uint32_t len = bam_cigar_oplen(cigar[i]);
            switch(op) {
                case 'N': // skipped region from the reference (intron)
                    enter_intron(op, len);
                    break;
                case '=':  // sequence match
                case 'M':  // alignment match (can be a sequence match or mismatch)
                case 'X':  // sequence mismatch
                    enter_aligned(op, len);
                    break;
                case 'D':  // deletion from the reference
                    enter_query_delete(op, len);
                    break;
                case 'I':  // insertion to the reference
                    enter_query_insert(op, len);
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                    enter_limbo();
                    break;
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    enter_limbo();
                    break;
                default:
                    throw std::invalid_argument("Unknown cigar operation '" + string(1, op) + "'");
            }
        }
        finish();
        *orientCnt = orientCnt_;
    }
};

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
    int first_strand = first_in_pair ^ reversed;
    int second_strand = second_in_pair ^ mate_reversed;
    if (strandness_ != UNSTRANDED) {
        // strandness_ is 0 for unstranded, 1 for RF, and 2 for FR
        int bool_strandness = strandness_ - 1;  // true for RF only
        first_strand ^= !bool_strandness;
        second_strand ^= !bool_strandness;
    }
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
    } else if (bam_aux_get(aln, "XS") != NULL) {
        return get_junction_strand_XS(aln);
    } else {
        // otherwise, we can't tell determine strand just from read orientation
        return '.';
    }
}


// Destructor
JunctionsExtractor::~JunctionsExtractor() {
    clear();
    if (out_sam_ != NULL) {
        sam_close(out_sam_);
    }
    bam_destroy1(aln_buf_);
    bam_hdr_destroy(in_header_);
    sam_close(in_sam_);
}


// determine if a junctions string in the for GT/AG is canonical.
// Do some basic qc on the junction
bool JunctionsExtractor::junction_qc(bam1_t *aln, hts_pos_t anchor_start, hts_pos_t intron_start,
                                     hts_pos_t intron_end, hts_pos_t anchor_end,
                                     uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt) {
    if ((aln->core).flag & BAM_FSECONDARY) {
        return false;
    } else if (((excludes_ & EXCLUDE_MULTI) != 0) && (get_num_aligns(aln) > 1)) {
        return false;
    } else if ((intron_end - intron_start < min_intron_length_) ||
               (intron_end - intron_start > max_intron_length_)) {
        return false;
    } else if (((intron_start - anchor_start) - left_mismatch_cnt < min_anchor_length_) ||
               ((anchor_end - intron_end) - right_mismatch_cnt < min_anchor_length_)) {
        return false;
    } else {
        return true;
    }
}

// Get the splice site strings, if sequence is not found in genome, return empty and
// output warning on first miss if instructed to ignore on error
string JunctionsExtractor::get_splice_sites(const string& chrom, char strand,
                                            hts_pos_t intron_start, hts_pos_t intron_end) {
    if (genome_ == NULL) {
        return "";
    }

    if (skip_missing_targets_ and not genome_->has_seq(chrom)) {
        if (missing_genomic_warned_.find(chrom) == missing_genomic_warned_.end()) {
            missing_genomic_warned_.insert(chrom);
            cerr << "Warning: genomic sequence not found for " << chrom << " splice junctions not available" << endl;
        }
        return "";
    }

    string splice_site = to_upper(genome_->fetch(chrom, intron_start, intron_start + 2) + "/"
                                  + genome_->fetch(chrom, intron_end - 2, intron_end));
    if (strand == '-') {
        splice_site = dna_reverse_complement(splice_site);
    }
    if (not is_canonical(splice_site)) {
        splice_site = to_lower(splice_site);
    }
    return splice_site;
}

Junction *JunctionsExtractor::create_junction(bam1_t *aln, const JunctionKey &key,
                                              const string& chrom, char strand,
                                              hts_pos_t anchor_start, hts_pos_t intron_start,
                                              hts_pos_t intron_end, hts_pos_t anchor_end) {
    char corrected_strand = strand;
    string splice_sites = get_splice_sites(chrom, strand, intron_start, intron_end);
    if ((splice_sites.size() > 0) && (strand == '.')) {
        // attempt to figure out strand from splice sites
        if (is_canonical(splice_sites)) {
            corrected_strand = '+';
        } else if (is_canonical(dna_reverse_complement(splice_sites))) {
            corrected_strand = '-';
            splice_sites = to_upper(dna_reverse_complement(splice_sites));
        }
    }

    if (trace_fh_ != NULL) {
        *trace_fh_ << chrom << "\t" << intron_start << "\t" << intron_end
                   << "\t" << bam_get_qname(aln) << "\t" << (aln->core).flag
                   << "\t" << strand << "\t" << corrected_strand
                   << "\t" << splice_sites
                   << "\t" << anchor_start << "\t" << anchor_end << endl;
    }
    Junction *junc = new Junction(chrom, intron_start, intron_end, anchor_start, anchor_end,
                                  corrected_strand, splice_sites, ijunc_);
    ijunc_++;
    junctions_[key] = junc;
    junc->count_read(get_read_category(aln), intron_start - aln->core.pos);
    return junc;
}

Junction *JunctionsExtractor::update_junction(bam1_t *aln, const JunctionKey &key,
                                              const string& chrom, char strand,
                                              hts_pos_t anchor_start, hts_pos_t intron_start,
                                              hts_pos_t intron_end, hts_pos_t anchor_end) {
    Junction *junc = junctions_[key];
    if (trace_fh_ != NULL) {
        *trace_fh_ << chrom << "\t" << intron_start << "\t" << intron_end
                   << "\t" << bam_get_qname(aln) << "\t" << (aln->core).flag
                   << "\t" << strand << "\t" << junc->strand
                   << "\t" << junc->splice_sites
                   << "\t" << anchor_start << "\t" << anchor_end << endl;
    }
    // Check if thick starts are any better
    if (anchor_start < junc->anchor_start)
        junc->anchor_start = anchor_start;
    if (anchor_end > junc->anchor_end)
        junc->anchor_end = anchor_end;
    junc->count_read(get_read_category(aln), intron_start - aln->core.pos);
    return junc;
}

// Add a junction to the junctions map
Junction *JunctionsExtractor::add_junction(bam1_t *aln, const string& chrom, char strand,
                                           hts_pos_t anchor_start, hts_pos_t intron_start,
                                           hts_pos_t intron_end, hts_pos_t anchor_end) {
    string splice_sites = get_splice_sites(chrom, strand, intron_start, intron_end);

    if (trace_fh_ != NULL) {
        *trace_fh_ << chrom << "\t" << intron_start << "\t" << intron_end
                   << "\t" << bam_get_qname(aln) << "\t" << (aln->core).flag
                   << "\t" << strand << "\t"
                   << "\t" << splice_sites
                   << "\t" << anchor_start << "\t" << anchor_end << endl;
    }
    JunctionKey key(chrom, intron_start, intron_end);

    // Check if new junction
    if (junctions_.count(key) == 0) {
        return create_junction(aln, key, chrom, strand, anchor_start, intron_start, intron_end, anchor_end);
    } else {
        return update_junction(aln, key, chrom, strand, anchor_start, intron_start, intron_end, anchor_end);
    }
}

// Validate a junction and save if it passes.
void JunctionsExtractor::record_junction(bam1_t *aln, const string& chrom, char strand,
                                         hts_pos_t anchor_start, hts_pos_t intron_start,
                                         hts_pos_t intron_end, hts_pos_t anchor_end,
                                         uint32_t left_mismatch_cnt, uint32_t right_mismatch_cnt,
                                         uint32_t left_indel_cnt, uint32_t right_indel_cnt,
                                         uint32_t left_indel_max, uint32_t right_indel_max,
                                         int *orientCnt) {
    assert((anchor_start >= 0) && (anchor_end >= 0));
    assert((intron_start >= 0) && (intron_end >= 0));
    assert(intron_start < intron_end);
    assert(anchor_start < anchor_end);
    assert((anchor_start < intron_start) && (intron_end < anchor_end));
    if (junction_qc(aln, anchor_start, intron_start, intron_end, anchor_end, left_mismatch_cnt, right_mismatch_cnt)) {
        Junction *junc = add_junction(aln, chrom, strand, anchor_start, intron_start, intron_end, anchor_end);
        if (junc->is_canonical()) {
            *orientCnt += (junc->strand == '+') ? 1 : -1;
        }
    }
}

// Process an alignment
void JunctionsExtractor::process_alignment(bam1_t *aln) {
    // skip if unmapped, only one cigar operation exists (likely all matches),
    // low-quality or duplicate
    if (!((aln->core.n_cigar <= 1) || (aln->core.tid < 0) || (aln->core.flag & (BAM_FQCFAIL | BAM_FDUP)))) {
        int orientCnt = 0;
        ReadJunctionExtractor readJuncExtract(aln, targets_[aln->core.tid], this);
        readJuncExtract.parse_read(&orientCnt);
        if (orientCnt != 0) {
            if (set_XS_strand_tag_) {
                update_strand_tag("XS", 'A', orientCnt, aln);
            }
            if (set_TS_strand_tag_) {
                update_strand_tag("TS", 'A', orientCnt, aln);
            }
        }
    }
}

// write aligment to pass through
void JunctionsExtractor::write_pass_through(bam1_t *aln,
                                            bam_hdr_t *in_header,
                                            samFile* out_sam) {
    if (out_sam != NULL) {
        if (sam_write1(out_sam, in_header, aln) < 0) {
            throw runtime_error("Error writing BAM record to pass-through file");
        }
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

// Open BAMs
void JunctionsExtractor::open(const string& bam,
                              const string& bam_pass_through) {
    bam_ = bam;
    in_sam_ = sam_open(bam_.c_str(), "r");
    if (in_sam_ == NULL) {
        throw runtime_error("Error opening BAM/SAM/CRAM file: " + bam_);
    }
    in_header_ = sam_hdr_read(in_sam_);
    save_targets(in_header_);
    if (trace_fh_ != NULL) {
        *trace_fh_ << "chrom" << "\t" << "intron_start" << "\t" << "intron_end"
                   << "\t" << "qname" << "\t" << "flag"
                   << "\t" << "strand" << "\t" << "corrected"
                   << "\t" << "splice_sites"
                   << "\t" << "anchor_start" << "\t" << "anchor_end" << endl;
    }
    if (bam_pass_through != "") {
        out_sam_ = open_pass_through(in_sam_, in_header_, bam_pass_through);
    }
    aln_buf_ = bam_init1();
}

// return pending or next
bam1_t* JunctionsExtractor::read_align() {
    if (aln_pending_) {
        aln_pending_ = false;
        return aln_buf_;
    } else {
        int stat = sam_read1(in_sam_, in_header_, aln_buf_);
        if (stat < -1) {
            throw runtime_error("Error reader BAM record: " + bam_);
        } else if (stat == -1) {
            return NULL;
        } else {
            return aln_buf_;
        }
    }
}

// check if an chrom has already been processed, indicating an unsorted BAM
void JunctionsExtractor::check_new_target(bam1_t *aln) {
    if (done_targets_.find(aln->core.tid) != done_targets_.end()) {
        throw runtime_error(targets_[aln->core.tid] + " already process, BAM appears to not be coordinate sorted");
    }
    done_targets_.insert(aln->core.tid);
}

// The workhorse - identifies junctions from BAM
void JunctionsExtractor::identify_junctions_for_target(int target_index) {
    bam1_t *aln;
    while((aln = read_align()) != NULL) {
        if (aln->core.tid != target_index) {
            aln_pending_ = true;
            check_new_target(aln);
            break;
        } else {
            process_alignment(aln);
            if (out_sam_ != NULL) {
                write_pass_through(aln, in_header_, out_sam_);
            }
        }
    }
}

// Copy remaining reads (unaligned) to pass-through, if it is open
void JunctionsExtractor::copy_unaligned_reads() {
    if (out_sam_ != NULL) {
        bam1_t *aln;
        while((aln = read_align()) != NULL) {
            write_pass_through(aln, in_header_, out_sam_);
        }
    }
}

