/* intronProspector.cc - command line program to extract 

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

#include <getopt.h>
#include <algorithm>
#include <fstream>
#include "junctions_extractor.hh"
#include "genome.hh"
#include "version.hh"

using namespace std;

static const uint32_t DEFAULT_MIN_ANCHOR_LENGTH = 8;
static const uint32_t DEFAULT_MIN_INTRON_LENGTH = 70;
static const uint32_t DEFAULT_MAX_INTRON_LENGTH = 500000;
static const bool DEFAULT_ALLOW_ANCHOR_INDELS = false;
static const uint32_t DEFAULT_MAX_ANCHOR_INDEL_SIZE = 36;
static const float DEFAULT_MIN_CONFIDENCE_SCORE = 0.0;
static const Strandness DEFAULT_STRANDED = UNSTRANDED;

static const char *usage_msg =
#include "intronProspector.man.h"
    ;

// Usage statement for this tool
static void usage() {
    cerr << usage_msg;
}

// convert a specification for strandness to constant.
static Strandness str_to_strandness(const string s) {
    // case-insensitive compare
    string su(s);
    transform(su.begin(), su.end(), su.begin(), ::toupper);
    if (su == "UN") {
        return UNSTRANDED;
    } else if (su == "RF") {
        return RF_STRANDED;
    } else if (su == "FR") {
        return FR_STRANDED;

    } else {
        cerr << "Error: invalid strandness value '" + s + "', expected one of 'UN', 'RF', or 'FR' (case insensitive)" << endl;
        exit(1);
    }
}

// parse an exclude category
static unsigned parse_exclude_cat(const string &cat) {
    if (cat == "multi") {
        return EXCLUDE_MULTI;
    } else {
        cerr << "Error: invalid exclude category '" << cat << "', expected one of: 'multi'" << endl;
        exit(1);
    }
}

// Parse command line
class CmdParser {
    public:
    // input
    string bam_file;
    uint32_t min_anchor_length;
    uint32_t min_intron_length;
    uint32_t max_intron_length;
    bool allow_anchor_indels;
    uint32_t max_anchor_indel_size;
    float min_confidence_score;
    Strandness strandness;
    unsigned excludes;
    string genome_fa;
    bool skip_missing_targets;

    // output
    string junction_bed;
    string intron_bed;
    string intron_call_tsv;
    string bam_pass_through;
    bool map_to_ucsc;
    string debug_trace_tsv;
    bool set_XS_strand_tag;
    bool set_TS_strand_tag;

    CmdParser(int argc, char *argv[]):
        bam_file("/dev/stdin"),
        min_anchor_length(DEFAULT_MIN_ANCHOR_LENGTH),
        min_intron_length(DEFAULT_MIN_INTRON_LENGTH),
        max_intron_length(DEFAULT_MAX_INTRON_LENGTH),
        allow_anchor_indels(DEFAULT_ALLOW_ANCHOR_INDELS),
        max_anchor_indel_size(DEFAULT_MAX_ANCHOR_INDEL_SIZE),
        min_confidence_score(DEFAULT_MIN_CONFIDENCE_SCORE),
        strandness(DEFAULT_STRANDED),
        excludes(EXCLUDE_NONE),
        skip_missing_targets(false),
        map_to_ucsc(false),
        set_XS_strand_tag(false),
        set_TS_strand_tag(false) {

        // definitions for long-only options
        static const int OPT_SET_XS_STRAND_TAG = 256;
        static const int OPT_SET_TS_STRAND_TAG = 257;
        
        struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"version", no_argument, NULL, 'v'},
            {"min-anchor-length", required_argument, NULL, 'a'},
            {"min-intron_length", required_argument, NULL, 'i'},
            {"max-intron_length", required_argument, NULL, 'I'},
            {"allow-anchor-indels", required_argument, NULL, 'd'},
            {"max-anchor-indel-size", required_argument, NULL, 'm'},
            {"min-confidence-score", required_argument, NULL, 'C'},
            {"strandness", required_argument, NULL, 's'},
            {"excludes", required_argument, NULL, 'X'},
            {"genome-fasta", required_argument, NULL, 'g'},
            {"skip-missing-targets", no_argument, NULL, 'S'},
            {"junction-bed", required_argument, NULL, 'j'},
            {"intron-bed", required_argument, NULL, 'n'},
            {"intron-calls", required_argument, NULL, 'c'},
            {"pass-through", required_argument, NULL, 'p'},
            {"map-to-ucsc", no_argument, NULL, 'U'},
            {"debug-trace", required_argument, NULL, 'D'},
            {"set-XS-strand-tag", no_argument, NULL, OPT_SET_XS_STRAND_TAG},
            {"set-TS-strand-tag", no_argument, NULL, OPT_SET_TS_STRAND_TAG},
            {NULL, 0, NULL, 0}
        };
            
        const char *short_options = "hva:i:I:C:s:X:g:S:j:n:c:p:UD:";
        int c;
        while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
            switch (c) {
                case 'a':
                    min_anchor_length = atoi(optarg);
                    break;
                case 'd':
                    allow_anchor_indels = true;
                    break;
                case 'm':
                    max_anchor_indel_size = atoi(optarg);
                    break;
                case 'i':
                    min_intron_length = atoi(optarg);
                    break;
                case 'I':
                    max_intron_length = atoi(optarg);
                    break;
                case 'C':
                    min_confidence_score = strtod(optarg, NULL);
                    break;
                case 's':
                    strandness = str_to_strandness(optarg);
                    break;
                case 'X':
                    excludes |= parse_exclude_cat(optarg);
                    break;
                case 'g':
                    genome_fa = string(optarg);
                    break;
                case 'S':
                    skip_missing_targets = true;
                    break;
                case 'j':
                    junction_bed = optarg;
                    break;
                case 'n':
                    intron_bed = optarg;
                    break;
                case 'c':
                    intron_call_tsv = optarg;
                    break;
                case 'p':
                    bam_pass_through = optarg;
                    break;
                case 'U':
                    map_to_ucsc = true;
                    break;
                case 'D':
                    debug_trace_tsv = optarg;
                    break;
                case OPT_SET_XS_STRAND_TAG:
                    set_XS_strand_tag = true;
                    break;
                case OPT_SET_TS_STRAND_TAG:
                    set_TS_strand_tag = true;
                    break;
                case 'h':
                    usage();
                    exit(0);
                case 'v':
                    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << " " << PACKAGE_URL << endl;
                    exit(0);
                case '?':
                default:
                    cerr << "Error: invalid option" << endl;
                    usage();
                    exit(1);
            }
        }

        if (argc - optind > 1) {
            cerr << "Error: too many positional arguments" << endl << endl;
            usage();
            exit(1);
        }
        if (argc - optind == 1) {
            bam_file = string(argv[optind++]);
        }
        if ((set_XS_strand_tag or set_TS_strand_tag)
            and ((bam_pass_through.size() == 0) or (genome_fa.size() == 0))) {
            cerr << "Error: --set-XS-strand-tag and --set-TS-strand-tag require --pass-through and --genome-fasta" << endl;
            exit(1);
        }
    }
};


static void output_junctions_for_target(JunctionsExtractor& extractor,
                                         float min_confidence_score,
                                         bool map_to_ucsc,
                                         ostream* junction_bed_fh,
                                         ostream* intron_bed_fh,
                                         ostream* intron_call_fh) {
    JunctionVector juncs = extractor.get_junctions();
    if (junction_bed_fh != NULL) {
        juncs.sort_by_anchors();
        print_anchor_bed(juncs, min_confidence_score, map_to_ucsc, *junction_bed_fh);
    }
    if ((intron_bed_fh != NULL) or (intron_call_fh != NULL)) {
        juncs.sort_by_introns();
        if (intron_bed_fh != NULL) {
            print_intron_bed(juncs, min_confidence_score, map_to_ucsc, *intron_bed_fh);
        }
        if (intron_call_fh != NULL) {
            print_intron_call_tsv(juncs, min_confidence_score, map_to_ucsc, *intron_call_fh);
        }
    }
}


static void extract_junctions(CmdParser &opts) {
    Genome *genome = NULL;
    if (opts.genome_fa.size() > 0) {
        genome = new Genome(opts.genome_fa);
    }
    ofstream *trace_fh = opts.debug_trace_tsv.size() > 0 ? new ofstream(opts.debug_trace_tsv.c_str()) :  NULL;
    JunctionsExtractor extractor(opts.min_anchor_length, opts.min_intron_length, opts.max_intron_length,
                                 opts.strandness, opts.excludes, genome,
                                 opts.skip_missing_targets,
                                 opts.set_XS_strand_tag, opts.set_TS_strand_tag,
                                 trace_fh);
    extractor.open(opts.bam_file, opts.bam_pass_through);
    ostream* junction_bed_fh = open_out_or_null(opts.junction_bed);
    ostream* intron_bed_fh = open_out_or_null(opts.intron_bed);
    ostream* intron_call_fh = open_out_or_null(opts.intron_call_tsv);
    if (intron_call_fh != NULL) {
        print_junction_call_header(*intron_call_fh);
    }
    for (int target_index = 0; target_index < extractor.get_num_targets(); target_index++) {
        extractor.identify_junctions_for_target(target_index);
        output_junctions_for_target(extractor, opts.min_confidence_score, opts.map_to_ucsc,
                                    junction_bed_fh, intron_bed_fh, intron_call_fh);
        extractor.clear();
    }
    extractor.copy_unaligned_reads();
    delete junction_bed_fh;
    delete intron_bed_fh;
    delete intron_call_fh;
    delete trace_fh;
    delete genome;
 }

// entry point
int main(int argc, char *argv[]) {
    CmdParser opts(argc, argv);
    try {
        extract_junctions(opts);
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << '\n';
        exit(1);
    }
    return 0;
}
