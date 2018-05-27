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
#include "junctions_extractor.hh"

using namespace std;

#if 0
//Usage statement for this tool
static int usage(const string& prog,
                 const string& msg) {
    
    out << "Usage:" 
        << "\t\t" << "" << endl;
    out << "Options:" << endl;
    out << "\t\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum \n"
        << "\t\t\t " << "anchor length on both sides are reported. [8]" << endl;
    out << "\t\t" << "-i INT\tMinimum intron length. [70]" << endl;
    out << "\t\t" << "-I INT\tMaximum intron length. [500000]" << endl;
    out << "\t\t" << "-o FILE\tThe file to write output to. [STDOUT]" << endl;
    out << "\t\t" << "-r STR\tThe region to identify junctions \n"
        << "\t\t\t " << "in \"chr:start-end\" format. Entire BAM by default." << endl;
    out << "\t\t" << "-s INT\tStrand specificity of RNA library preparation \n"
        << "\t\t\t " << "(0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). [1]" << endl;
    out << endl;
    return 0;
}

// Parse the options passed to this tool
int main(int argc, char *argv[]) {
    optind = 1; // Reset before parsing again.
    int c;
    min_anchor_length = 8;
    min_intron_length = 70;
    max_intron_length = 500000;
    string bam_file = "/dev/stdin";
    string output_file = "/dev/stdout";
    while((c = getopt(argc, argv, "ha:i:I:o:r:s:")) != -1) {
        switch(c) {
            case 'a':
                min_anchor_length = atoi(optarg);
                break;
            case 'i':
                min_intron_length = atoi(optarg);
                break;
            case 'I':
                max_intron_length = atoi(optarg);
                break;
            case 'o':
                output_file = string(optarg);
                break;
#if 0 // FIXME
            case 'r':
                region_ = string(optarg);
                break;
#endif
            case 's':
                strandness_ = atoi(optarg);
                break;
            case 'h':
                usage("");
            case '?':
            default:
                usage();
                exit(1)
        }
    }
    if (argc - optind >= 1) {
        bam_ = string(argv[optind++]);
    }
    if (optind < argc || bam_ == "NA") {
        usage();
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    return 0;
}
#else
// FIXME: quick hack
int main(int argc, char *argv[]) {
    JunctionsExtractor je("../tests/input/test1.sam");
    je.identify_junctions_from_bam();
    je.print_all_junctions(cout);
}

#endif
