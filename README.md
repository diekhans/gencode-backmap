# gencode-backmap
Mapping of GENCODE gene annotation set files to older assembies

This provides tools to map the genome corrdinates in various GENCODE download files to the previous release of the reference
genome.  This uses UCSC liftOver change and the base-by-base transmap algorithm to produce mapping.  Reports are made on genes
that don't map exactly.
