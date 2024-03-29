# gencode-backmap
Package for mapping of GENCODE gene annotation files to older assemblies.

This provides tools to map the genome coordinates in GENCODE to previous
releases of the reference genome.  This projects annotations through genomics
alignment using the TransMap algorithm to produce mappings.  Reports are made
on status of all mappings.

This maps all GENCODE GTF and GFF3 files. It does not map the PolyA_feature
file. It is recommended to only use these on primary chromosomes, due to the
high-similarity of alt-loci sequences to the primary assembly coupled with the
growth in number of alt-loci making correct back-mapping difficult.

## Description

### Algorithm

The program takes the current (source) GENCODE GFF3 or GTF, cross-assembly
genomic alignments, and the previous (target) GENCODE annotations.

Mapping is done on a per-gene mapping using the following steps:

- Project transcripts of the gene through the alignments, keeping
  exons chained.
  - If there are multiple mappings, first look for ones that the overlapping
    to the previous version of the transcript, if it exists.  Otherwise, if
    there is a previous version of the gene, select mappings overlapping the
    gene.  Otherwise, to filter for paralog mappings, pick the mapping with
    the most similar span as the source.
  - Project features of the transcript, such as CDS and start codons,
    the transcript alignment between the genomes.  This ensures that
    features stay in the same location within the transcript.
- Check all transcripts of the gene for consistency.  Reject source gene mappings
  with transcripts on different chromosomes or strands, or where
  the genomic length of the gene has changed more than 50%.
- If a version of the gene exists in the target and the mapped gene doesn't
  overlap the target gene, it is also rejected.
  - If a gene did not map or was rejected and a version of the gene with the
    same biotype exists in the target annotations, use the existing gene.
- Small, automatic-only or all automatic genes are optionally not mapped, with
  the target annotation being passed through This avoids complex mappings of
  small RNAs imported from other database (e.g. mirRNAs).
- Target genes with no corresponding mappings and that overlap patched regions
  or regions with GRC indecent reports in the target genome may optionally be
  passed through.  This address a fair number of problem cases.  This was a
  common problem on GRCh37 chrX.

Pairing of source and target genes is somewhat complex due to instability of
some gene identifiers between assemblies.  If a matching base gene id (less
version) is not found, an attempt is made to match the genes using the
symbolic name.

### Identification

Gene and transcript identifiers (`*_id` attribute) are based on the source
identifiers.  For annotations that are passed through unchanged from the
target assembly, the identifiers are using as-is.  For mapped annotations, the
a mapping version number, in the form `_n`, is appended to the id.  Where `n`
is one based number indicating a version of the mapping.  The mapping version
start at `1' and are addition to the standard identifier version, for example
`ENST00000456328.2_1`

When a previous release file is supplied, the gene and transcript annotations
are checked and the mapping version incremented if defined attributes change.

Genes have the mapping version numbers incremented if any of the following
change from the previous release:
- contained transcripts, include version and mapping version number
- coordinates of the gene change
- the gene_type or gene_status changes
- the gene_name changes

Transcripts have the mapping version numbers incremented if any of the
following change from the previous release:
- coordinates of the transcript or it's features change
- the transcript_type or transcript_status changes
- the transcript_name change
- any tag changes
- CDS or other sub-features of a transcript change.

The following ids have to the mapping versions applied
- `gene_id`
- `transcript_id`
- 'havana_gene
- 'havana_transcript

To deal with the, the PAR used the name stable ids (until V44), mapping are
restricted to the same chromosome.

### Categorization of mappings 

Information is collected on mappings and saved as attributes of the
GFF3/GTF records.  They are also records at the gene and transcript level
the mapping information file.  The attributes and their values are:

- `remap_status` - Attribute that indicates the status of the mapping. Possible values are:
  - `full_contig` - Gene or transcript completely mapped to the target genome with all features intact.
  - `full_fragment` -  Gene or transcript completely to the target genome with insertions in some features.
      These are usually small insertions.
  - `partial` - Gene or transcript partially mapped to the target genome.
  - `deleted` - Gene or transcript did not map to the target genome.
  - `no_seq_map` - The source sequence is not in the assembly alignments. This
      will occur with alt loci genes if the alignments only contain the primary assembly.
  - `gene_conflict` - Transcripts in the gene mapped to multiple locations.
  - `gene_size_change` - Transcripts caused gene's length to change by more than 50%.
     This is to detect mapping to processed pseudogenes and mapping across tandem gene duplications.
  - `automatic_small_ncrna_gene` - Gene is a from a small, automatic (ENSEMBL source) non-coding RNA.  These
     are take from the target annotations if `--useTargetForAutoSmallNonCoding` is specified.
  - `automatic_gene` - Gene is a from an automatic process (ENSEMBL source).  These
     are take from the target annotations if `--useTargetForAutoGenes` is specified.
     It is not recommended to use --useTargetForAutoGenes, instead use
     --useTargetForAutoSmallNonCoding to exclude the automatic small ncRNAs, to
     get the other automatic genes mapped.
  - `pseudogene` - Pseuduogene annotations (excluding polymorphic).  These
     are take from the target annotations if `--useTargetForPseudoGenes` is specified.
     It is not recommended to use --useTargetForPseudoGenes.
- `remap_original_id` - Original ID attribute of the feature.  If a feature is split when mapped,
  new IDs are created, otherwise the original ID is used.
- `remap_original_location` - Location of the feature in the source genome.
- `remap_num_mappings` - Number of mappings of the feature (on the same chromosome), only one of them was used.
- `remap_target_status` - Attribute that compares the mapping to the existing target annotations. Possible values are:
  - `new` - Gene or transcript was not in target annotations.
  - `lost` -Gene or transcript exists in source and target genome, however source was not mapped.
  - `overlap` - Gene or transcript overlaps previous version of annotation on target genome.
  - `nonOverlap` - Gene or transcript exists in target, however source mapping is to a different location.
    This is often mappings to a gene family members or pseudogenes.
- `remap_substituted_missing_target` - target gene annotate was substituted

### Requirements
- a reasonable modern C++ compiler, tested with g++ (GCC) 4.8.5
- python3, recommended running from a virtual environment
  ```
  python3 -m virtualenv gencode-venv
  source gencode-venvvenv/bin/activate
  pip install -r  requirements.txt
  ```

### Usage

The following files are needed to map using the UCSC liftover alignments:

- UCSC liftover alignments.  These are the alignments of GRCh37 to GRCh38,
  which will be swapped.  This produced slightly better results that using the
  GRCh38 to GRCh37 alignments.
  - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
- NCBI assembly report file for GRCh37
  - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt
- NCBI assembly report file for GRCh38.  Since only the primary assembly is
  map, it is not necessary to update this for new releases.
  - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_report.txt
- GENCODE V19 primary annotations in GFF3 and GTF
  - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
  - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
- Current GENCODE primary assembly GFF3 and GTF files to map
  - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz
  - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.basic.annotation.gtf.gz
- GRC patch BED file for GRCh37 from UCSC.  Download from
  - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hg19Patch13Patches.txt.gz
- GRC incident table from UCSC browser:
  - http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/grcIncidentDb/hg19.grcIncidentDb.bb
- Previous GFF3 from (e.g. gencode.v28lift37.annotation.gff3.gz )


Remove bin column from patches and convert BigBed to bed and combine.
This does not need to be redone each releases.

```
   zcat hg19Patch13Patches.txt.gz | cut -f 2- > hg19Patch13Patches.bed
   bigBedToBed Hg19.grcIncidentDb.bb Hg19.grcIncidentDb.bed
   cat hg19Patch13Patches.bed Hg19.grcIncidentDb.bed >problemRegions.bed
```

GENCODE uses different sequence naming for unplaced chromosome sequences than
UCSC.  Additionally, GENCODE annotates the GRCh37-lite chrM, while UCSC had
incorporated a different chrM version.  Convert the UCSC alignments to
use GENCODE names:
```
../gencode-backmap/bin/ucscLiftEdit hg38ToHg19.over.chain.gz GCF_000001405.36_GRCh38.p10_assembly_report.txt GCF_000001405.13_GRCh37_assembly_report.txt hg38ToHg19.over.gencode.chain
```

Map the annotation files:
```
../gencode-backmap/bin/gencode-backmap --swapMap --useTargetForAutoSmallNonCoding --onlyManualForTargetSubstituteOverlap --substituteMissingTargets=V19 --headerFile=liftGxfHeader.txt --targetGxf=gencode.v19.annotation.gff3.gz --targetPatches=problemRegions.bed  --previousMappedGxf=gencode.v28lift37.annotation.gff3.gz  gencode.v29.annotation.gff3.gz hg38ToHg19.over.gencode.chain gencode.v29lift37.annotation.gff3 gencode.v29lift37.map-info.tsv
../gencode-backmap/bin/gencode-backmap --swapMap --useTargetForAutoSmallNonCoding --onlyManualForTargetSubstituteOverlap --substituteMissingTargets=V19 --headerFile=liftGxfHeader.txt --targetGxf=gencode.v19.annotation.gtf.gz  --targetPatches=problemRegions.bed  --previousMappedGxf=gencode.v28lift37.annotation.gff3.gz gencode.v29.annotation.gtf.gz   hg38ToHg19.over.gencode.chain gencode.v29lift37.annotation.gtf
```

Where `liftGxfHeader.txt` is the comments to add at the beginning of the output GFF3 or GTF files.
This does not include GFF3 meta comment.

### Installation

#### Requirements
- GCC supporting -std=c++11 (4.7 or newer), or clang
- UCSC Browser library (kent/src)
- to run the tests, the following programs from the UCSC Browser are required:
  - gff3ToGenePred
  - gtfToGenePred
- Python 3.7

### Compiling

- Obtain and compile the libraries and command line utilities for the kent library:
  - http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob;f=src/product/README.building.source

- Edit `config.mk` to set the following variables
  - `KENTDIR` - point to kent/src directory
  - `SAMTABIXDIR` - if browser library is compiled with samtabix support, this
    is used to find the library.
- Compile code with `make` and turn tests with `make test`
- There is no install step, use directly from the bin directory


### Version numbering
- The file `src/version.h.in` contains the version number.  It should be
  update when doing a new release.  The git sha1 hash of the commit is also added
  to this file. Be sure to do a make clean before doing a final build to get
  this file rebuild before a release.

- A tag is created in the form v1.0_rel on master.  Due to the limited number
  of users of the code, release branches are not created.


