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
#include <vector>
#include "junctions_extractor.hh"
#include "strconvert.hh"
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
    exit(0);
}

// error and reference to -h
static void cmd_error(const string& msg) {
    cerr << "Error: "  << msg << endl;
    cerr << "Use --help to get usage" << endl;
    exit(1);
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
    vector<string> bam_files;
    uint32_t min_anchor_length;
    uint32_t min_intron_length;
    uint32_t max_intron_length;
    bool allow_anchor_indels;
    bool unsorted;
    uint32_t max_anchor_indel_size;
    float min_confidence_score;
    Strandness strandness;
    unsigned excludes;
    string genome_fa;
    bool skip_missing_targets;

    // output
    string junction_bed;
    string intron_bed;
    string intron_bed6;
    string intron_call_tsv;
    string bam_pass_through;
    string debug_trace_tsv;
    bool set_XS_strand_tag;
    bool set_TS_strand_tag;

    CmdParser(int argc, char *argv[]):
        min_anchor_length(DEFAULT_MIN_ANCHOR_LENGTH),
        min_intron_length(DEFAULT_MIN_INTRON_LENGTH),
        max_intron_length(DEFAULT_MAX_INTRON_LENGTH),
        allow_anchor_indels(DEFAULT_ALLOW_ANCHOR_INDELS),
        max_anchor_indel_size(DEFAULT_MAX_ANCHOR_INDEL_SIZE),
        min_confidence_score(DEFAULT_MIN_CONFIDENCE_SCORE),
        strandness(DEFAULT_STRANDED),
        excludes(EXCLUDE_NONE),
        skip_missing_targets(false),
        set_XS_strand_tag(false),
        set_TS_strand_tag(false) {

        try {
            parse_cmd_args(argc, argv);
        } catch (const std::exception& ex) {
            cerr << "Error parsing command line: " << ex.what() << endl;
            exit(1);
        }
    }


    private:
    void parse_cmd_args(int argc, char *argv[]) {
        // definitions for long-only options
        static const int OPT_SET_XS_STRAND_TAG = 256;
        static const int OPT_SET_TS_STRAND_TAG = 257;
        
        struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"version", no_argument, NULL, 'v'},
            {"unsorted", no_argument, NULL, 'u'},
            {"min-anchor-length", required_argument, NULL, 'a'},
            {"min-intron-length", required_argument, NULL, 'i'},
            {"max-intron-length", required_argument, NULL, 'I'},
            {"allow-anchor-indels", no_argument, NULL, 'd'},
            {"max-anchor-indel-size", required_argument, NULL, 'm'},
            {"min-confidence-score", required_argument, NULL, 'C'},
            {"strandness", required_argument, NULL, 's'},
            {"excludes", required_argument, NULL, 'X'},
            {"genome-fasta", required_argument, NULL, 'g'},
            {"skip-missing-targets", no_argument, NULL, 'S'},
            {"junction-bed", required_argument, NULL, 'j'},
            {"intron-bed", required_argument, NULL, 'n'},
            {"intron-bed6", required_argument, NULL, 'b'},
            {"intron-calls", required_argument, NULL, 'c'},
            {"pass-through", required_argument, NULL, 'p'},
            {"debug-trace", required_argument, NULL, 'D'},
            {"set-XS-strand-tag", no_argument, NULL, OPT_SET_XS_STRAND_TAG},
            {"set-TS-strand-tag", no_argument, NULL, OPT_SET_TS_STRAND_TAG},
            {NULL, 0, NULL, 0}
        };
            
        const char *short_options = "hvua:i:I:C:s:X:g:S:j:n:b:c:p:UD:";
        int c;
        while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
            switch (c) {
                case 'a':
                    min_anchor_length = atoi(optarg);
                    break;
                case 'u':
                    unsorted = true;
                    break;
                case 'd':
                    allow_anchor_indels = true;
                    break;
                case 'm':
                    max_anchor_indel_size = toUnsigned(optarg);
                    break;
                case 'i':
                    min_intron_length = toUnsigned(optarg);
                    break;
                case 'I':
                    max_intron_length = toUnsigned(optarg);
                    break;
                case 'C':
                    min_confidence_score = toFloat(optarg);
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
                case 'b':
                    intron_bed6 = optarg;
                    break;
                case 'c':
                    intron_call_tsv = optarg;
                    break;
                case 'p':
                    bam_pass_through = optarg;
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
                case 'v':
                    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << " " << PACKAGE_URL << endl;
                    exit(0);
                case '?':
                default:
                    cmd_error("invalid option");
            }
        }

        // if indels are not allowed, also indicate this my setting max size
        if (not allow_anchor_indels) {
            max_anchor_indel_size = 0;
        }

        while (optind < argc) {
            bam_files.push_back(string(argv[optind++]));
        }
        if (bam_files.size() == 0) {
            bam_files.push_back("/dev/stdin");
        }
        if ((bam_files.size() > 1) and (bam_pass_through != "")) {
            cmd_error("can only have one input BAM with --bam-pass-through");
        }
        if ((set_XS_strand_tag or set_TS_strand_tag)
            and ((bam_pass_through.size() == 0) or (genome_fa.size() == 0))) {
            cmd_error("--set-XS-strand-tag and --set-TS-strand-tag require --pass-through and --genome-fasta");
        }
    }
};


class OutputFiles {
    public:
    const string bam_pass_through;
    ostream* junction_bed_fh;
    ostream* intron_bed_fh;
    ostream* intron_bed6_fh;
    ostream* intron_call_fh;

    OutputFiles(const CmdParser& opts):
        bam_pass_through(opts.bam_pass_through),
        junction_bed_fh(open_out_or_null(opts.junction_bed)),
        intron_bed_fh(open_out_or_null(opts.intron_bed)),
        intron_bed6_fh(open_out_or_null(opts.intron_bed6)),
        intron_call_fh(open_out_or_null(opts.intron_call_tsv)) {
        if (intron_call_fh != NULL) {
            print_junction_call_header(*intron_call_fh);
        }
    }

    ~OutputFiles() {
        delete junction_bed_fh;
        delete intron_bed_fh;
        delete intron_bed6_fh;
        delete intron_call_fh;
    }
};

static void output_junctions(JunctionsExtractor& extractor,
                             float min_confidence_score,
                             OutputFiles& output) {
    JunctionVector juncs = extractor.get_junctions();
    if (output.junction_bed_fh != NULL) {
        juncs.sort_by_anchors();
        print_anchor_bed(juncs, min_confidence_score, *output.junction_bed_fh);
    }
    if ((output.intron_bed_fh != NULL) or (output.intron_bed6_fh != NULL)
        or (output.intron_call_fh != NULL)) {
        juncs.sort_by_introns();
        if (output.intron_bed_fh != NULL) {
            print_intron_bed(juncs, min_confidence_score, 9, *output.intron_bed_fh);
        }
        if (output.intron_bed6_fh != NULL) {
            print_intron_bed(juncs, min_confidence_score, 6, *output.intron_bed6_fh);
        }
        if (output.intron_call_fh != NULL) {
            print_intron_call_tsv(juncs, min_confidence_score, *output.intron_call_fh);
        }
    }
}

/* Extract junctions for the single BAM, optimizing memory,
 * by process one target at a time.
 */
static void extract_junctions_by_target(JunctionsExtractor& extractor,
                                        const string& bam_file,
                                        float min_confidence_score,
                                        OutputFiles& output) {
    extractor.open(bam_file);
    if (output.bam_pass_through != "") {
        extractor.open_pass_through(output.bam_pass_through);
    }
    for (int target_index = 0; target_index < extractor.get_num_targets(); target_index++) {
        extractor.identify_junctions_for_target(target_index);
        output_junctions(extractor, min_confidence_score, output);
        extractor.clear();
    }
    extractor.copy_unaligned_reads();
    extractor.close();
}

/* Extract junctions for an unsorted BAM.  This allows pass-through.
 */
static void extract_junctions_unsorted(JunctionsExtractor& extractor,
                                       const string& bam_file,
                                       float min_confidence_score,
                                       OutputFiles& output) {
    extractor.open(bam_file);
    if (output.bam_pass_through != "") {
        extractor.open_pass_through(output.bam_pass_through);
    }
    extractor.identify_junctions_for_bam();
    extractor.copy_unaligned_reads();
    extractor.close();
    output_junctions(extractor, min_confidence_score, output);
}
                                
/* Extract junctions for the multiple BAM files.
 */
static void extract_junctions_multiple(JunctionsExtractor& extractor,
                                       const vector<string>& bam_files,
                                       float min_confidence_score,
                                       OutputFiles& output) {
    for (int ibam = 0; ibam < bam_files.size(); ibam++) {
        extractor.open(bam_files[ibam]);
        extractor.identify_junctions_for_bam();
        extractor.close();
    }
    output_junctions(extractor, min_confidence_score, output);
}
                                
static void extract_junctions(CmdParser &opts) {
    Genome *genome = NULL;
    if (opts.genome_fa.size() > 0) {
        genome = new Genome(opts.genome_fa);
    }
    ofstream *trace_fh = opts.debug_trace_tsv.size() > 0 ? new ofstream(opts.debug_trace_tsv.c_str()) :  NULL;
    JunctionsExtractor extractor(opts.min_anchor_length, opts.max_anchor_indel_size,
                                 opts.min_intron_length, opts.max_intron_length,
                                 opts.strandness, opts.excludes, genome,
                                 opts.skip_missing_targets,
                                 opts.set_XS_strand_tag, opts.set_TS_strand_tag,
                                 trace_fh);
    OutputFiles output(opts);
    if (opts.bam_files.size() > 1) {
        extract_junctions_multiple(extractor, opts.bam_files,
                                   opts.min_confidence_score, output);
    } else if (opts.unsorted) {
        extract_junctions_unsorted(extractor, opts.bam_files[0],
                                   opts.min_confidence_score, output);
    } else {
        extract_junctions_by_target(extractor, opts.bam_files[0],
                                    opts.min_confidence_score, output);
    }
    delete trace_fh;
    delete genome;
 }

// entry point
int main(int argc, char *argv[]) {
    CmdParser opts(argc, argv);
    try {
        extract_junctions(opts);
    } catch (const std::exception& e) {
        cmd_error(e.what());
    }
    return 0;
}
