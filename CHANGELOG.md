# Major and user-visible changes

## version 1.4.0 [in progress]
* Add --sj-filter=canon option `intronProspector` to keep only canonical splice junctions.
  This is the default if a genome is provide.  Filtering can be disabled with
  --sj-filter=all.
* Add --sj-filter=canon option `intronProspectorMerge` to keep only canonical splice
  junctions.


## version 1.3.0 2025-04-10
* Added ability to read and write gzip compressed TSV and BED files based on an name
  extension of '.gz'

## version 1.2.0 2025-03-15
* Added ability to process multiple BAMs at once.  This will use more memory
  and possibly improve the confidence scoring in some cases
* Added option to allow input BAM to be unsorted.  This will use more memory.
