/*  genome.cc -- splice junction extraction from genome.

    Copyright (c) 2019, Mark Diekhans, University of California, Santa Cruz

**/
#include "genome.hh"
#include "type_ops.hh"
#include <stdexcept>
#include <stdlib.h>
#include "htslib/faidx.h"

/* reverse DNA, index by base */
static char *reverseTbl = NULL;

static void buildReverseTbl() {
    reverseTbl = new char[128];
    // fill in with pass-through first
    for (int i = 0; i < 128; i++) {
        reverseTbl[i] = i;
    }
    reverseTbl['A'] = 'T';
    reverseTbl['T'] = 'A';
    reverseTbl['G'] = 'C';
    reverseTbl['C'] = 'G';
    reverseTbl['a'] = 't';
    reverseTbl['t'] = 'a';
    reverseTbl['g'] = 'c';
    reverseTbl['c'] = 'g';
}

/* reverse-complement a DNA string */
string dnaReverseComplement(const string& dna) {
    if (reverseTbl == NULL) {
        buildReverseTbl();
    }
    string rcDna;
    rcDna.resize(dna.size());
    int j = dna.size() - 1;
    for (int i = 0; i < dna.size(); i++) {
        rcDna[j--] = reverseTbl[dna[i]];
    }
    return rcDna;
}

Genome::Genome(const string& genome_fa):
    faidx(NULL) {
    faidx = fai_load3(genome_fa.c_str(), NULL, NULL, 0);
    if (faidx == NULL) {
        throw runtime_error("can't open genome FASTA file or associated index: " + genome_fa);
    }
}

Genome::~Genome() {
    if (faidx != NULL) {
        fai_destroy(faidx);
    }
}

bool Genome::has_seq(const string &chrom) const {
    return faidx_has_seq(faidx, chrom.c_str());
}

const string Genome::fetch(const string &chrom, int start, int end) const {
    int retlen = 0;
    char *seq = faidx_fetch_seq(faidx, chrom.c_str(), start, end - 1, &retlen);
    if (seq == NULL) {
        throw runtime_error("can't load genome sequence for: " + mk_coords_str(chrom, start, end));
    }
    if (retlen != (end - start)) {
        throw runtime_error("unexpected length returned for: "  + mk_coords_str(chrom, start, end) + ": " + int_to_string(retlen));
    }
    string sseq(seq);
    free(seq);
    return sseq;
}

