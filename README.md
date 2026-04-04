# intron-prospector
Identify putative introns in RNA-Seq alignments in BAM/SAM alignments.
This works with both short and long-read RNA sequencing.

- `intron-prospector` - calls introns from SAM/BAM files
- `intron-prospector-merge` - merges output from multiple `intron-prospector` runs or convert between formats

Usage is described in the [intron-prospector manual page](docs/intron-prospector.md)
and [intron-prospector-merge manual page](docs/intron-prospector-merge.md)


## Installation

[See INSTALL.md for installation](INSTALL.md)

The only dependencies are `htslib` and a modern C++ compiler

## Origins

This code is derived from the Griffith Lab's excellent 
[`RegTools`](https://github.com/griffithlab/regtools) package.
