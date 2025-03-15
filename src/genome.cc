/*  genome.cc -- splice junction extraction from genome.

    Copyright (c) 2019-2025, Mark Diekhans, University of California, Santa Cruz

**/
#include "genome.hh"
#include "type_ops.hh"
#include <stdexcept>
#include <stdlib.h>

/* reverse DNA, index by base */
static char *reverse_tbl = NULL;

static void build_reverse_tbl() {
    reverse_tbl = new char[128];
    // fill in with pass-through first
    for (int i = 0; i < 128; i++) {
        reverse_tbl[i] = i;
    }
    reverse_tbl[int('A')] = 'T';
    reverse_tbl[int('T')] = 'A';
    reverse_tbl[int('G')] = 'C';
    reverse_tbl[int('C')] = 'G';
    reverse_tbl[int('a')] = 't';
    reverse_tbl[int('t')] = 'a';
    reverse_tbl[int('g')] = 'c';
    reverse_tbl[int('c')] = 'g';
}

/* reverse-complement a DNA string */
string dna_reverse_complement(const string& dna) {
    if (reverse_tbl == NULL) {
        build_reverse_tbl();
    }
    string rcDna;
    rcDna.resize(dna.size());
    int j = dna.size() - 1;
    for (int i = 0; i < dna.size(); i++) {
        rcDna[j--] = reverse_tbl[int(dna[i])];
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

