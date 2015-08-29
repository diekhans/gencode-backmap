# gencode-backmap
Mapping of GENCODE gene annotation set files to older assembies

This provides tools to map the genome corrdinates in various GENCODE download
files to the previous release of the reference genome.  It uses UCSC liftOver
alignments and the base-by-base transmap algorithm to produce mappings.
Reports are made on genes that don't map exactly or are dropped.

The maps all GENCODE GTF and GFF3 files. It does not map the PolyA_feature
file.

### Usage

~/compbio/code/pycbio-package/bin/ncbiAssemblyReportConvert --fromId=refSeqAccn sizes data/GCF_000001405.25.assembly.txt data/GRCh37.p13.sizes
ftp://ftp.ncbi.nlm.nih.gov/pub/murphyte/assm_alignments/GRCh38.p2-GRCh37.p13.gff.gz

ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.assembly.txt
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.assembly.txt

bin/ncbiAssemblyReportConvert --fromId=refSeqAccn sizes data/GCF_000001405.25.assembly.txt data/GRCh37.p13.sizes
bin/ncbiAssemblyReportConvert --fromId=refSeqAccn sizes data/GCF_000001405.28.assembly.txt data/GRCh38.p2.sizes



### Installation

#### Requirements
- Python 2.7
- GCC supporting -std=c++11 (4.7 or newer), or clang

