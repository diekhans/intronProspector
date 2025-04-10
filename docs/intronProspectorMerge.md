# NAME

**intronProspectorMerge** — Merge putative introns junctions calls made `intronProspector` and/or convert output formats

# SYNOPSIS

`intronProspectorMerge [options] [intronCallsTsv ...]`

# DESCRIPTION

Merge output of the `intronProspector` intron calls tab-separated (TSV) files, as created by the `--intron-calls` option.
This program can also be used to convert the output format from a single run.

Compressed intron calls TSV files are recognized if they end in `.gz`.
TSV and BED files will be automatically compressed with `gzip` if they end in `.gz`.

# Options

`-h, --help`

> Prints brief usage information.

`-v, --version`

> Prints the current version number.

`-i FILE, --intron-calls-files=FILE`

> File containing a list input intron calls files, one file per line'

`-c FILE, --intron-calls=FILE`

> Write a tab-separated values (TSV) file, with the junctions calls.  The confidence column will be the maximum confidence in any of the files.

`-j FILE, --junction-bed=FILE`

> Write the junction calls and support anchors to this file.  This is in the same format as TopHat `junctions.bed` and `regtools junction extract` output.  It is a UCSC BED track, with each junction consists of two connected BED blocks, where each block is as long as the maximal overhang of any read spanning the junction. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

`-n FILE, --intron-bed=FILE`

> Write the intron BED with the bounds of the intron. The score is the number of alignments spanning the junction, with a maximum score of 1000 for UCSC browser compatibility.  If genome is supplied, BED is colored green for U2 junctions, blue for U12, or red for unknown.

`-b FILE, --intron-bed6=FILE`

> Write the intron BED 6 with the bounds of the introns. The score is the number of alignments spanning the junction.   This is for software not wanting to create a browser track.

# BUGS

See GitHub Issues: <https://github.com/diekhans/intronProspector/issues>

# AUTHOR

Mark Diekhans <markd@ucsc.edu>

Source available from <https://github.com/diekhans/intronProspector>

Base on code from RegTools <https://github.com/griffithlab/regtools>
by Avinash Ramu <aramu@genome.wustl.edu>.

