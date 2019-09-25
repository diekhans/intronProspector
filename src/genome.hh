/*  genome.hh -- splice junction extraction from genome.

    Copyright (c) 2019, Mark Diekhans, University of California, Santa Cruz

**/
#ifndef GENOME_HH
#define GENOME_HH
#include <string>

typedef struct __faidx_t faidx_t;

using namespace std;

/** class used to load and cache splice junctions from genome FASTA */
class Genome {
    private:
    faidx_t *faidx;
    public:
    Genome(const string& genomeFa);
    ~Genome();

    const string fetch(const string &chrom, int start, int end);
};

/* reverse-complement a DNA string */
string dnaReverseComplement(const string& dna);

#endif
