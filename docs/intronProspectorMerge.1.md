NAME
====

**intronProspectorMerge** — Merge putative introns junctions calls made `intronProspector`

SYNOPSIS
========

`intronProspectorMerge [options] intronCallsTsv ...`

DESCRIPTION
===========

Merge output of the `intronProspector` intron calls tab-separated (TSV) files, as created by the `--intron-calls` option.


Options
-------

`-h, --help`

> Prints brief usage information.

`-v, --version`

> Prints the current version number.

`-U, --map-to-ucsc`

> Naively generate UCSC chromosome names in TSV and BED files.  This pre-pends `chr` to numeric and X/Y names and change `MT` to `chrMT`, other name are not modified.  This will not produce the correct results for other sequences such as alts, patches, and unmapped sequences.  It does not modify records passed through (-p).

`-c FILE, --intron-calls=FILE`

> Write a tab-separated values (TSV) file, with the junctions calls.  The confidence column will be the maximum confidence in any of the files.

`-j FILE, --junction-bed=FILE`

> Write the junction calls and support anchors to this file.  This is in the same format as TopHat `junctions.bed` and `regtools junction extract` output.  It is a UCSC BED track, with each junction consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

`-n FILE, --intron-bed=FILE`

> Write the intron BED with the bounds of the intron. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

BUGS
====

See GitHub Issues: <https://github.com/diekhans/intronProspector/issues>

AUTHOR
======

Mark Diekhans <markd@ucsc.edu>

Source available from <https://github.com/diekhans/intronProspector>

Base on code from RegTools <https://github.com/griffithlab/regtools>
by Avinash Ramu <aramu@genome.wustl.edu>.

