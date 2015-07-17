# gencode-backmap
Mapping of GENCODE gene annotation set files to older assembies

This provides tools to map the genome corrdinates in various GENCODE download
files to the previous release of the reference genome.  It uses UCSC liftOver
alignments and the base-by-base transmap algorithm to produce mappings.
Reports are made on genes that don't map exactly or are dropped.

The maps all GENCODE GTF and GFF3 files, as well as the PolyA_feature files.

### Installation



