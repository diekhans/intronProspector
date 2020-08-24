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
