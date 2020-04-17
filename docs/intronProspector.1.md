NAME
====

**intronProspector** — Extract putative junctions from RNA-Seq alignments

SYNOPSIS
========

`intronProspector [options] [readaligns]`

DESCRIPTION
===========

Find putative intron junctions in a RNA-Seq alignment. The *readaligns* file maybe in SAM, BAM, or CRAM format and does not need to be sorted or indexed and maybe streamed. If omitted, stdin is used.

This program allows for calling splice junctions indebtedness of the alignment program.  It maybe used in a pipeline, copying the alignment file on `stdin` to `stdout`.  It can sit between an aligner outputting SAM and `samtools` to convert to BAM/CRAM.

Options
-------

`-h, --help`

> Prints brief usage information.

`-v, --version`

> Prints the current version number.

`-a INT, --min-anchor-length=INT`

> Minimum anchor length. Junctions which satisfy a minimum anchor length on both sides are reported.  Mismatch bases don't count towards the meeting this threshold.  The default is 8 bases.

`-i INT, --min-intron_length=INT`

> Minimum intron length. The default is 70 bases.

`-I INT,  --max-intron_length=INT`

> Maximum intron length. The default is 500000 bases.

`-C FLOAT, --min-confidence-score=FLOAT`

> Calculate the Shannon-Wiener Diversity Index to use as a confidence score and discard intron calls below this value.  The default is 0.0, which discards no calls. A value of 1.0 is a good threshold for filtering. This methodology it taken from JuncBASE DOI: 10.1101/gr.108662.110.

`-s STRING, --strandness=STRING`

> Strand specificity of RNA library preparation.  Use `UN` for unstranded, `RF` for first-strand, `FR` for second-strand (case-insensitive).  The default is `UN`.  This is used to set the strand in the junction-format BED file.

`-g fasta, --genome-fasta=fasta`

> Genome FASTA file, must be indexed by `samtools faidx`.  Donor and acceptor dinuclotides are identified if provided. 

`-S, --skip-missing-targets`

> Skip getting splice junctions when target sequence is missing in genome FASTA rather than generate an error and stop.

`-U, --map-to-ucsc`

> Naively generate UCSC chromosome names in TSV and BED files.  This pre-pends `chr` to numeric and X/Y names and change `MT` to `chrMT`, other name are not modified.  This will not produce the correct results for other sequences such as alts, patches, and unmapped sequences.  It does not modify records passed through (-p).

`-c FILE, --intron-calls=FILE`

> Write a tab-separated values (TSV) file, with the junctions calls.  It will contain the following columns :
* chrom - chromosome name
* intron_start - zero based starting coordinates of intron.
* intron_end - zero based ending coordinates of intron.
* strand - `+`, `-`, or `.` if not known.  This is set only if experiment is stranded, it is not set from genomic sequence.
* uniq_mapped_count - number of uniquely mapped reads supporting the junction.
* multi_mapped_count - number of uniquely mapped reads supporting the junction.
* unsure_mapped_count - number of reads that are either in discordant or partial mapped pair or do not have the BAM `NH` tag.
* max_left_overhang - maximum number of bases overlapping the exon upstream of the intron.
* max_right_overhang - maximum number of bases overlapping the exon downstream of the intron.
* confidence - confidence score.
* splice_sites - Donor-acceptor pair, if genome provided. In upper-case for known splicing patterns, lower-case for unknown.  If strand is not known, this field is empty.

`-j FILE, --junction-bed=FILE`

> Write the junction calls and support anchors to this file.  This is in the same format as ToHat `junctions.bed` and `regtools junction extract` output.  It UCSC BED track, with each junction consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.

`-n FILE, --intron-bed=FILE`

> Write the intron BED with the bounds of the intron. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.

`-p FILE, --pass-through=FILE`

> Pass through input BAM/SAM records to this file, used for constructing pipelines with `/dev/stdout` is specified.  CRAM output is not support, these will be written as BAM.

`-D FILE, --debug-trace=FILE`

> Output records, in TSV format, for reach read intron indicating the information going into classifying it, including read name.  First few columns are BED-like for easy conversion.

`-X category, --exclude=category`

> Exclude reads or introns from this category.  Current categories are:
> *multi* - excluding multi-mapped reads.

NOTES
=====
The computation of strand is problematic.  If the strandness of the experiment is specified, then that is used to determine stand.  If the alignment provides the XS attribute, that is used.  Otherwise, the strand can't be determined from the BAM.  If the genome is provided and a known splice sites are detected, this is then used if the stand is not identified by other methods.

Secondary alignments are not used to support introns.

EXAMPLES
========

Call junctions from a BAM file, also creating BEDs of junctions and introns:
```
intronProspector --intron-calls=introns.tsv --junction-bed=juncs.bed --intron-bed=introns.bed reads.bam
```

Pipeline to call introns and create a CRAM file:
```
cat reads.sam \
    | samtools sort -O sam  \
    | ./intronProspector -c introns.tsv -p /dev/stdout \
    | samtools view -O CRAM -T grch38.fa >reads.cram
```
Note that the `cat` command could be an aligner output a SAM file and that the genome FASTA file must be index by `samtools faidx`.


BUGS
====

See GitHub Issues: <https://github.com/diekhans/intronProspector/issues>

AUTHOR
======

Mark Diekhans <markd@ucsc.edu>

Source available from <https://github.com/diekhans/intronProspector>

Base on code from RegTools <https://github.com/griffithlab/regtools>
by Avinash Ramu <aramu@genome.wustl.edu>.

