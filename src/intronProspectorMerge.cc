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
#include <string>
#include <vector>
#include "junctions.hh"
#include "tsv.hh"
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

    // output
    string junction_bed;
    string intron_bed;
    string intron_call_tsv;

    CmdParser(int argc, char *argv[]) {
        
        struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"version", no_argument, NULL, 'v'},
            {"junction-bed", required_argument, NULL, 'j'},
            {"intron-bed", required_argument, NULL, 'n'},
            {"intron-calls", required_argument, NULL, 'c'},
            {NULL, 0, NULL, 0}
        };
            
        const char *short_options = "hvj:n:c:U";
        int c;
        while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
            switch (c) {
                case 'j':
                    junction_bed = optarg;
                    break;
                case 'n':
                    intron_bed = optarg;
                    break;
                case 'c':
                    intron_call_tsv = optarg;
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

        if (argc - optind < 1) {
            cmd_error("a least one input call TSV required");
        }
        while (optind < argc) {
            input_calls_tsvs.push_back(string(argv[optind++]));
        }
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

static void process_junction(JunctionTable& junction_tbl,
                             const Junction& junc) {
    JunctionKey key(junc.chrom, junc.intron_start, junc.intron_end);
    if (junction_tbl.count(key) == 0) {
        junction_tbl[key] = new Junction(junc);
    } else {
        junction_tbl[key]->merge(junc);
    }
}


static void process_calls_tsv(JunctionTable& junction_tbl,
                              const string& tsv_file) {
    Tsv tsv(tsv_file);
    while (tsv.next_row()) {
        process_junction(junction_tbl, read_junction(tsv));
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
        process_calls_tsv(junction_tbl, opts.input_calls_tsvs[i]);
    }

    // output
    JunctionVector juncs = junction_tbl.get_junctions();
    juncs.sort_by_anchors();
    renumber_junctions(juncs);
    if (opts.junction_bed != "") {
        ofstream out(opts.junction_bed.c_str());
        print_anchor_bed(juncs, 0.0, out);
    }
    if ((opts.intron_bed != "") or (opts.intron_call_tsv != "")) {
        juncs.sort_by_introns();
        if (opts.intron_bed != "") {
            ofstream out(opts.intron_bed.c_str());
            print_intron_bed(juncs, 0.0, out);
        }
        if (opts.intron_call_tsv != "") {
            ofstream out(opts.intron_call_tsv.c_str());
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
