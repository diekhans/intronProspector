/*  genome.hh -- splice junction extraction from genome.

    Copyright (c) 2019-2025, Mark Diekhans, University of California, Santa Cruz

**/
#ifndef GENOME_HH
#define GENOME_HH
#include <string>

#include "htslib/faidx.h"

using namespace std;

/** class used to load and cache splice junctions from genome FASTA */
class Genome {
    private:
    faidx_t *faidx;
    public:
    Genome(const string& genomeFa);
    ~Genome();

    bool has_seq(const string &chrom) const;

    const string fetch(const string &chrom, int start, int end) const;
};

/* reverse-complement a DNA string */
string dna_reverse_complement(const string& dna);

#endif
