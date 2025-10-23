# intronProspector
Identify putative introns in RNA-Seq alignments in BAM/SAM alignments.
This works with both short and long-read RNA sequencing.

- `intronProspector` - calls introns from SAM/BAM files
- `intronProspectorMerge` - merges output from multiple `intronProspector` runs or convert between formats

Usage is described in the [intronProspector manual page](docs/intronProspector.md)
and [intronProspectorMerge manual page](docs/intronProspectorMerge.md)


## Installation

[See INSTALL.md for instation](INSTALL.md)

The only dependencies are `htslib` and a modern C++ compiler

## Origins

This code is derived from the Griffith Lab's excellent 
[`RegTools`](https://github.com/griffithlab/regtools) package.
