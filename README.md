# gencode-backmap
Mapping of GENCODE gene annotation set files to older assembies

This provides tools to map the genome corrdinates in various GENCODE download
files to the previous release of the reference genome.
alignments and the base-by-base TransMap algorithm to produce mappings.
Reports are made on genes that don't map exactly or are dropped.

The maps all GENCODE GTF and GFF3 files. It does not map the PolyA_feature
file.

### Usage

  cd data/
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/murphyte/assm_alignments/GRCh38.p2-GRCh37.p13.gff.gz
  wget -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.assembly.txt
  wget -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.assembly.txt
  cd ..
  ../bin/ncbiAssemblyReportConvert --fromId=refSeqAccn sizes data/GCF_000001405.25.assembly.txt data/GRCh37.p13.sizes
  ../bin/ncbiAssemblyReportConvert --fromId=refSeqAccn sizes data/GCF_000001405.28.assembly.txt data/GRCh38.p2.sizes

../bin/ncbiAssemblyReportConvert --fromIdType=refSeqAccn --toIdType=gencode lift data/GCF_000001405.25.assembly.txt data/GRCh37.p13.refseq-gencode.lift
../bin/ncbiAssemblyReportConvert --fromIdType=refSeqAccn --toIdType=gencode lift data/GCF_000001405.28.assembly.txt data/GRCh38.p2.refseq-gencode.lift

 query, main GFF3 sequence columns, is GRCh37, target, in attributes, is GRCh38
 gff3ToPsl data/GRCh37.p13.sizes data/GRCh38.p2.sizes data/GRCh38.p2-GRCh37.p13.gff.gz data/GRCh38.p2-GRCh37.p13.refseq.psl
  liftUp -nohead -type=.psl /dev/stdout data/GRCh38.p2.refseq-gencode.lift error data/GRCh38.p2-GRCh37.p13.refseq.psl | liftUp -nohead -type=.psl -pslQ data/GRCh38.p2-GRCh37.p13.gencode.psl data/GRCh37.p13.refseq-gencode.lift error /dev/stdin

### Output
- tags

### Installation

#### Requirements
- Python 2.7
- GCC supporting -std=c++11 (4.7 or newer), or clang
- kent library
 gff3ToGenePred
 gtfToGenePred
 gff3ToPsl 


### Alignments
http://www.ncbi.nlm.nih.gov/genome/tools/remap/docs/alignments
ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

For assembly-assembly alignments:
- reciprocity=3 == alignment is best for both Reference and Target. Also referred to as 'First-Pass Alignments'
- reciprocity=1 or 2 == alignment is best for one sequence, but a better alignment is reported for the other sequence. Also referred to as 'Second-Pass Alignments'
