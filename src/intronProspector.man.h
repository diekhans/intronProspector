"\n"
"\n"
"NAME\n"
"\n"
"\n"
"INTRONPROSPECTOR — Extract putative intron junctions from RNA-Seq\n"
"alignments\n"
"\n"
"\n"
"\n"
"SYNOPSIS\n"
"\n"
"\n"
"intronProspector [options] [readaligns]\n"
"\n"
"\n"
"\n"
"DESCRIPTION\n"
"\n"
"\n"
"Find putative intron junctions in a RNA-Seq alignment. The _readaligns_\n"
"file maybe in SAM, BAM, or CRAM format and does not need to be sorted or\n"
"indexed and maybe streamed. If omitted, stdin is used.\n"
"\n"
"This program allows for integrating splice junction calling into an\n"
"alignment pipeline. Pass-through mode copys the alignment file on stdin\n"
"to stdout. It can then sit between an aligner outputting SAM and\n"
"samtools that converts to BAM/CRAM.\n"
"\n"
"Both short-read (Illumina) and long-read RNA-Seq can be process. For\n"
"short-reads, it is recommended to use --min-confidence-score=1.0.\n"
"Introns are as determined by the aligner and indicated in the BAM by N\n"
"operations.\n"
"\n"
"\n"
"Options\n"
"\n"
"-h, --help\n"
"\n"
"  Prints brief usage information.\n"
"\n"
"-v, --version\n"
"\n"
"  Prints the current version number.\n"
"\n"
"-a INT, --min-anchor-length=INT\n"
"\n"
"  Minimum anchor length. Junctions which satisfy a minimum anchor length\n"
"  on both sides are reported. Mismatch bases don’t count towards the\n"
"  meeting this threshold. The default is 8 bases.\n"
"\n"
"-i INT, --min-intron_length=INT\n"
"\n"
"  Minimum intron length. The default is 70 bases.\n"
"\n"
"-I INT,  --max-intron_length=INT\n"
"\n"
"  Maximum intron length. The default is 500000 bases.\n"
"\n"
"-d,  --allow-anchor-indels\n"
"\n"
"  Are indels allowed in anchors? Useful for ONT reads. Indels bases\n"
"  don’t count towards the meeting the anchor size threshold. The default\n"
"  is don’t allow indels. .\n"
"\n"
"-m INT,  --max-anchor-indel-size=INT\n"
"\n"
"  Maximum size of any contiguous indel in an anchor if indels are\n"
"  allowed in anchors. Exceeding this length discards the intron. The\n"
"  default is 36.\n"
"\n"
"-C FLOAT, --min-confidence-score=FLOAT\n"
"\n"
"  Calculate the Shannon-Wiener Diversity Index to use as a confidence\n"
"  score and discard intron calls below this value. The default is 0.0,\n"
"  which discards no calls. A value of 1.0 is a good threshold for\n"
"  filtering. This methodology it taken from JuncBASE DOI:\n"
"  10.1101/gr.108662.110.\n"
"\n"
"-s STRING, --strandness=STRING\n"
"\n"
"  Strand specificity of RNA library preparation. Use UN for unstranded,\n"
"  RF for first-strand, FR for second-strand (case-insensitive). The\n"
"  default is UN. This is used to set the strand in the junction-format\n"
"  BED file.\n"
"\n"
"-g fasta, --genome-fasta=fasta\n"
"\n"
"  Genome FASTA file, must be indexed by samtools faidx. Donor and\n"
"  acceptor dinuclotides are identified if provided.\n"
"\n"
"-S, --skip-missing-targets\n"
"\n"
"  Skip getting splice junctions when target sequence is missing in\n"
"  genome FASTA rather than generate an error and stop.\n"
"\n"
"-c FILE, --intron-calls=FILE\n"
"\n"
"  Write a tab-separated values (TSV) file, with the junctions calls. It\n"
"  will contain the following columns : * chrom - chromosome name *\n"
"  intron_start - zero based starting coordinates of intron. * intron_end\n"
"  - zero based ending coordinates of intron. * strand - +, -, or . if\n"
"  not known. This is set only if experiment is stranded, it is not set\n"
"  from genomic sequence. * uniq_mapped_count - number of uniquely mapped\n"
"  reads supporting the junction. * multi_mapped_count - number of\n"
"  uniquely mapped reads supporting the junction. * unsure_mapped_count -\n"
"  number of reads that are either in discordant or partial mapped pair\n"
"  or do not have the BAM NH tag. * max_left_overhang - maximum number of\n"
"  bases overlapping the exon upstream of the intron. *\n"
"  max_right_overhang - maximum number of bases overlapping the exon\n"
"  downstream of the intron. * confidence - confidence score. *\n"
"  splice_sites - Donor-acceptor pair, if genome provided. In upper-case\n"
"  for known splicing patterns, lower-case for unknown. If strand is not\n"
"  known, this field is empty.\n"
"\n"
"-j FILE, --junction-bed=FILE\n"
"\n"
"  Write the junction calls and support anchors to this BED 12 file. This\n"
"  is in the same format as TopHat junctions.bed and\n"
"  regtools junction extract output. It is a UCSC BED track, with each\n"
"  junction consists of two connected BED blocks, where each block is as\n"
"  long as the maximal overhang of any read spanning the junction. The\n"
"  score is the number of alignments spanning the junction, with a\n"
"  maximum score of 1000 for UCSC browser compatibility. If genome is\n"
"  supplied, BED is colored green for U2 junctions, blue for U12, or red\n"
"  for unknown.\n"
"\n"
"-n FILE, --intron-bed=FILE\n"
"\n"
"  Write the intron BED 9 with the bounds of the introns. The score is\n"
"  the number of alignments spanning the junction, with a maximum score\n"
"  of 1000 for UCSC browser compatibility. If genome is supplied, BED is\n"
"  colored green for U2 junctions, blue for U12, or red for unknown.\n"
"\n"
"-p FILE, --pass-through=FILE\n"
"\n"
"  Pass through input BAM/SAM records to this file, used for constructing\n"
"  pipelines with /dev/stdout is specified. CRAM output is not support,\n"
"  these will be written as BAM.\n"
"\n"
"-D FILE, --debug-trace=FILE\n"
"\n"
"  Output records, in TSV format, for reach read intron indicating the\n"
"  information going into classifying it, including read name. First few\n"
"  columns are BED-like for easy conversion.\n"
"\n"
"-X category, --exclude=category\n"
"\n"
"  Exclude reads or introns from this category. Current categories are:\n"
"  _multi_ - excluding multi-mapped reads.\n"
"\n"
"--set-XS-strand-tag\n"
"\n"
"  Set the XS:A tag to the strand if stand for a read can be determined\n"
"  from the introns based on a count of recognized vs unrecognized splice\n"
"  sites. If the strand can’t be determined, the tags are not modified,\n"
"  possible leaving an existing XS:A tag in place. Requires\n"
"  --pass-through and --genome-fasta.\n"
"\n"
"--set-TS-strand-tag\n"
"\n"
"  Set the TS:A tag in the same manner as --set-XS-strand-tag.\n"
"\n"
"\n"
"\n"
"NOTES\n"
"\n"
"\n"
"The computation of strand is problematic. If the strandness of the\n"
"experiment is specified, then that is used to determine stand. If the\n"
"alignment provides the XS attribute, that is used. Otherwise, the strand\n"
"can’t be determined from the BAM. If the genome is provided and a known\n"
"splice sites are detected, this is then used if the stand is not\n"
"identified by other methods.\n"
"\n"
"Secondary alignments are not used to support introns.\n"
"\n"
"\n"
"\n"
"EXAMPLES\n"
"\n"
"\n"
"Call junctions from a BAM file, also creating BEDs of junctions and\n"
"introns:\n"
"\n"
"    intronProspector --intron-calls=introns.tsv --junction-bed=juncs.bed --intron-bed=introns.bed reads.bam\n"
"\n"
"Pipeline to call introns and create a CRAM file:\n"
"\n"
"    cat reads.sam \\n"
"        | samtools sort -O sam  \\n"
"        | ./intronProspector -c introns.tsv -p /dev/stdout \\n"
"        | samtools view -O CRAM -T grch38.fa >reads.cram\n"
"\n"
"Note that the cat command could be an aligner output a SAM file and that\n"
"the genome FASTA file must be index by samtools faidx.\n"
"\n"
"\n"
"\n"
"BUGS\n"
"\n"
"\n"
"See GitHub Issues: https://github.com/diekhans/intronProspector/issues\n"
"\n"
"\n"
"\n"
"AUTHOR\n"
"\n"
"\n"
"Mark Diekhans markd@ucsc.edu\n"
"\n"
"Source available from https://github.com/diekhans/intronProspector\n"
"\n"
"Base on code from RegTools https://github.com/griffithlab/regtools by\n"
"Avinash Ramu aramu@genome.wustl.edu.\n"
