/* intronProspectorMerge.cc - merge from multiple intronProspector runs

    Copyright (c) 2018-2025, Mark Diekhans, University of California, Santa Cruz

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
#include <string>
#include <vector>
#include "zfstream.hh"
#include "junctions.hh"
#include "tsv.hh"
#include "splice_juncs.hh"
#include "version.hh"

using namespace std;

static const char *usage_msg =
#include "intronProspectorMerge.man.h"
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

// Parse command line
class CmdParser {
    public:
    // input
    vector<string> input_calls_tsvs;
    string intron_calls_files;

    // output
    string junction_bed;
    string intron_bed;
    string intron_bed6;
    string intron_call_tsv;

    // filter
    JunctionFilter junction_filter;
    
    CmdParser(int argc, char *argv[]):
        junction_filter(NULL_SJ_FILTER) {
        try {
            parse_cmd_args(argc, argv);
        } catch (const std::exception& ex) {
            cerr << "Error parsing command line: " << ex.what() << endl;
            exit(1);
        }
    }

    private:
    int parse_options(int argc, char *argv[]) {
        struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"version", no_argument, NULL, 'v'},
            {"intron-calls-files", required_argument, NULL, 'i'},
            {"junction-bed", required_argument, NULL, 'j'},
            {"intron-bed", required_argument, NULL, 'n'},
            {"intron-bed6", required_argument, NULL, 'b'},
            {"intron-calls", required_argument, NULL, 'c'},
            {"sj-filter", required_argument, NULL, 'f'},
            {NULL, 0, NULL, 0}
        };
            
        const char *short_options = "hvi:j:n:b:c:f:";
        int c;
        while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
            switch (c) {
                case 'i':
                    intron_calls_files = optarg;
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
                case 'f':
                    junction_filter = junction_filter_parse(string(optarg));
                    break;
                case 'h':
                    usage();
                case 'v':
                    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << " " << PACKAGE_URL << endl;
                    exit(0);
                case '?':
                default:
                    cmd_error("bug: option not handled");
            }
        }
        return optind;
    }

    void read_input_call_list(const string& intron_calls_files) {
        AutoGzipInput infh(intron_calls_files);
        string line;
        while (getline(infh, line)) {
            input_calls_tsvs.push_back(line);
        }
    }

    void collect_input_calls(int argc, char *argv[], int nextopt) {
        while (nextopt < argc) {
            input_calls_tsvs.push_back(string(argv[nextopt++]));
        }
        if (intron_calls_files.size() > 0) {
            read_input_call_list(intron_calls_files);
        }
        if (input_calls_tsvs.size() == 0) {
            cmd_error("a least one input call TSV required on the command list or with the --intron_call_files (-i) option");
        }
    }

    void parse_cmd_args(int argc, char *argv[]) {
        int nextopt = parse_options(argc, argv);
        collect_input_calls(argc, argv, nextopt);
    }
};

// current row as Junction
static Junction read_junction(const Tsv& tsv) {
    // ijunc is ignored
    uint32_t start = tsv.get_col_int("intron_start");
    uint32_t end = tsv.get_col_int("intron_end");
    Junction junc(tsv.get_col("chrom"), start, end,
                  start - tsv.get_col_int("max_left_overhang"),
                  end + tsv.get_col_int("max_right_overhang"),
                  tsv.get_col("strand")[0], tsv.get_col("splice_sites"), 0);
    junc.set_read_counts(SINGLE_MAPPED_READ, tsv.get_col_int("uniq_mapped_count"));
    junc.set_read_counts(MULTI_MAPPED_READ, tsv.get_col_int("multi_mapped_count"));
    junc.set_read_counts(UNSURE_READ, tsv.get_col_int("unsure_mapped_count"));
    junc.set_confidence(tsv.get_col_float("confidence"));
    return junc;
}

static void process_junction(JunctionFilter junction_filter,
                             JunctionTable& junction_tbl,
                             const Junction& junc) {
    JunctionKey key(junc.chrom, junc.intron_start, junc.intron_end);
    if (junction_tbl.count(key) == 0) {
        junction_tbl[key] = new Junction(junc);
    } else {
        junction_tbl[key]->merge(junc);
    }
}

static void process_calls_tsv(JunctionFilter junction_filter,
                              JunctionTable& junction_tbl,
                              const string& tsv_file) {
    Tsv tsv(tsv_file);
    while (tsv.next_row()) {
        if (junction_filter_check(junction_filter, tsv.get_col("splice_sites"))) {
            process_junction(junction_filter, junction_tbl, read_junction(tsv));
        }
    }
}

static void renumber_junctions(JunctionVector& juncs) {
    for (int i = 0; i < juncs.size(); i++) {
        juncs[i]->ijunc = i;
    }
}

static void intron_prospector_merge(CmdParser &opts) {
    // load and merge
    JunctionTable junction_tbl;
    for (int i = 0; i < opts.input_calls_tsvs.size(); i++) {
        process_calls_tsv(opts.junction_filter, junction_tbl, opts.input_calls_tsvs[i]);
    }

    // output
    JunctionVector juncs = junction_tbl.get_junctions();
    juncs.sort_by_anchors();
    renumber_junctions(juncs);
    if (opts.junction_bed != "") {
        AutoGzipOutput out(opts.junction_bed);
        print_anchor_bed(juncs, 0.0, out);
    }
    if ((opts.intron_bed != "") or (opts.intron_bed6 != "") or (opts.intron_call_tsv != "")) {
        juncs.sort_by_introns();
        if (opts.intron_bed != "") {
            AutoGzipOutput out(opts.intron_bed);
            print_intron_bed(juncs, 0.0, 9, out);
        }
        if (opts.intron_bed6 != "") {
            AutoGzipOutput out(opts.intron_bed6);
            print_intron_bed(juncs, 0.0, 6, out);
        }
        if (opts.intron_call_tsv != "") {
            AutoGzipOutput out(opts.intron_call_tsv);
            print_junction_call_header(out);
            print_intron_call_tsv(juncs, 0.0, out);
        }
    }
 }

// entry point
int main(int argc, char *argv[]) {
    CmdParser opts(argc, argv);
    try {
        intron_prospector_merge(opts);
    } catch (const std::exception& e) {
        cmd_error(e.what());
    }
    return 0;
}
