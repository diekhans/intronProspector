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

#include "junctions.hh"
#include <fstream>
#include "type_ops.hh"
#include <assert.h>
#include <stdexcept>

using namespace std;

float Junction::NULL_CONFIDENCE = nanf("not-a-number");

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

string Junction::get_description() const {
    return mk_coords_str(chrom, intron_start, intron_end) + "[" + strand + "/" + splice_sites + "]";
}

void Junction::merge(const Junction& j1) {
    if ((j1.chrom != chrom) || (j1.intron_start != intron_start) || (j1.intron_end != intron_end)) {
        throw logic_error("attempt to merge different introns:" + j1.get_description() + " with " + get_description());
    }
    if (j1.strand != strand) {
        throw runtime_error("intron records have different strand:" + j1.get_description() + " with " + get_description());
    }
    if (j1.splice_sites != splice_sites) {
        throw runtime_error("intron records have different splice sites:" + j1.get_description() + " with " + get_description());
    }
    anchor_start = min(j1.anchor_start, anchor_start);
    anchor_end = max(j1.anchor_end, anchor_end);
    read_counts[SINGLE_MAPPED_READ] += j1.read_counts[SINGLE_MAPPED_READ];
    read_counts[MULTI_MAPPED_READ] += j1.read_counts[MULTI_MAPPED_READ];
    read_counts[UNSURE_READ] += j1.read_counts[UNSURE_READ];
    confidence = max(j1.confidence, confidence);
}

// Compare two junctions
bool JunctionVector::junctions_lt(const Junction *j1,
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

// Compare two introns
bool JunctionVector::introns_lt(const Junction *j1,
                                const Junction *j2) {
    if (j1->chrom < j2->chrom){
        return true;
    }
    if (j1->chrom > j2->chrom){
        return false;
    }
    // Same chromosome
    if (j1->intron_start < j2->intron_start) {
        return true;
    }
    if (j1->intron_start > j2->intron_start) {
        return false;
    }
    // Same start
    if (j1->intron_end < j2->intron_end) {
        return true;
    }
    if (j1->intron_end > j2->intron_end) {
        return false;
    }
    // Same end
    return j1->strand < j2->strand;
}

// Trim read counts to fit in BED score restriction or 0..1000
static unsigned read_count_to_bed_score(uint64_t read_count) {
    return (read_count <= 1000) ? read_count : 1000;
}

// naive mapping to ucsc chrom name
static const string make_ucsc_chrom(const string& chrom) {
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
static const string& getBedColor(const string& splice_sites) {
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
void Junction::print_anchor_bed(unsigned ijunc,
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
        << "\t" << getBedColor(splice_sites) << "\t" << 2
        << "\t" << intron_start - anchor_start << "," << anchor_end - intron_end
        << "\t" << "0," << intron_end - anchor_start << endl;
}

// Print BED with intron as block  ijunc is used to
// make the BED name.
void Junction::print_intron_bed(unsigned ijunc,
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
        << "\t" << intron_start << "\t" << intron_end << "\t" << getBedColor(splice_sites) << endl;
}

// Print header for junction call TSV
void Junction::print_junction_call_header(ostream& out) {
    out << "chrom" << "\t" << "intron_start" << "\t" << "intron_end" << "\t" << "strand"
        << "\t" << "uniq_mapped_count" << "\t" << "multi_mapped_count" << "\t" << "unsure_mapped_count"
        << "\t" << "max_left_overhang" << "\t" << "max_right_overhang" << "\t" << "confidence"
        << "\t" << "splice_sites" << endl;
}

// Print row to junction call TSV
void Junction::print_junction_call_row(bool map_to_ucsc,
                                       ostream& out) const {
    if (map_to_ucsc) {
        out << make_ucsc_chrom(chrom);
    } else {
        out << chrom;
    }        
    out << "\t" << intron_start << "\t" << intron_end << "\t" << strand
        << "\t" << read_counts[SINGLE_MAPPED_READ] << "\t" << read_counts[MULTI_MAPPED_READ] << "\t" << read_counts[UNSURE_READ]
        << "\t" << intron_start - anchor_start << "\t" << anchor_end - intron_end
        << "\t" << get_confidence() << "\t" << splice_sites << endl;
}

// Print BED with anchors as blocks and intron as gap.
void print_anchor_bed(const JunctionVector& juncs,
                      float min_confidence_score,
                      bool map_to_ucsc,
                      const string& outfile) {
    ofstream out(outfile.c_str());
    for (unsigned ijunc = 0; ijunc < juncs.size(); ijunc++) {
        if (juncs[ijunc]->get_confidence() >= min_confidence_score) {
            juncs[ijunc]->print_anchor_bed(ijunc, map_to_ucsc, out);
        }
    }
}

// Print BED with intron as block
void print_intron_bed(const JunctionVector& juncs,
                      float min_confidence_score,
                      bool map_to_ucsc,
                      const string& outfile) {
    ofstream out(outfile.c_str());
    for (unsigned ijunc = 0; ijunc < juncs.size(); ijunc++) {
        if (juncs[ijunc]->get_confidence() >= min_confidence_score) {
            juncs[ijunc]->print_intron_bed(ijunc, map_to_ucsc, out);
        }
    }
}

// Print TSV with intron information
void print_intron_call_tsv(const JunctionVector& juncs,
                           float min_confidence_score,
                           bool map_to_ucsc,
                           const string& outfile) {
    ofstream out(outfile.c_str());
    Junction::print_junction_call_header(out);
    out.precision(3);
    for (unsigned ijunc = 0; ijunc < juncs.size(); ijunc++) {
        if (juncs[ijunc]->get_confidence() >= min_confidence_score) {
            juncs[ijunc]->print_junction_call_row(map_to_ucsc, out);
        }
    }
}

