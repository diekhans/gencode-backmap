ROOT = ..
include ${ROOT}/config.mk

### IMPORTANT: test GFF3, if new are added changes thests
GRCh37_GFF3s = data/gencode.v19.annotation.gff3 data/idChange.v19.gff3 \
	data/dup-transcripts.v19.gff3
GRCh38_GFF3s = data/gencode.v22.annotation.gff3 data/idChange.v25.gff3 data/gencode.v29.annotation.gff3 \
	data/gencode.v31.annotation.gff3 data/dup-transcripts.v33.gff3


targetGtfArg = --targetGxf=data/gencode.v19.annotation.gtf
targetGff3Arg = --targetGxf=data/gencode.v19.annotation.gff3
targetSubstArg = --substituteMissingTargets=V19
targetPatchArg = --targetPatches=data/problemRegions.bed
headerArg = --headerFile=data/header.txt

valgrind = valgrind --leak-check=full --num-callers=50

testUcscLiftOverChains = data/hg38ToHg19.over.chain
testGencodeLiftOverChains = output/hg38ToHg19.gencode.over.chain
testNcbiLiftOverChains = data/GRCh38.p2-GRCh37.p13.gencode.psl

diff = diff -u

# use transcript_id attr, not ID
gff3ToGenePred = gff3ToGenePred -geneNameAttr=gene_id  -rnaNameAttr=transcript_id

# commands to normalize genePreds for comparison because the V19 GTF's from
# Ensembl at least sometimes doens't include a start_codon so whack the CDS
# completeness columns. Also filter genePred to remove tests sequences were
# gtfToGenePred joins adjacent blocks and gff3ToGenePred doesn't
# also wack PAR transcripts because they get a differnt name with GTF.
genePredCmpIgnore = ENST00000578049|ENST00000624585
genePredParAdjust = gsub("ENSTR", "ENST0", $$0); gsub("ENSGR", "ENSG0", $$0); gsub("_PAR_Y", "", $$0)
normalizeGenePred =  awk 'BEGIN{IF=OFS="\t"} /${genePredCmpIgnore}/{next} {${genePredParAdjust}; $$13=$$14="foo"; print $$0}' | sort
ncbiAssemblyReportConvert = ../bin/ncbiAssemblyReportConvert

# compare by genePred ensure that GTF and GFF3 produce the same results

all: test

test: gff3UcscTest gtfUcscTest cmpUcscTest \
	gff3ParNamingTest gtfParNamingTest cmpParNamingTest \
	gff3NcbiTest gtfNcbiTest \
	gff3UcscSubstituteAuto gff3UcscSubstituteAutoSmallNcRna \
	gff3UcscSubstituteManOverlap gtfUcscSubstituteManOverlap \
	cmpUcscSubstituteManOverlap \
	gff3V29RegressTest \
	reportsTests ucscLiftEditTest regressTests

gff3UcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gff3 data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${gff3ToGenePred} output/$@.unmapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.unmapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfUcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gtf data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.unmapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.unmapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.unmapped.gtf output/$@.unmapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

cmpUcscTest: gff3UcscTest gtfUcscTest
	${diff} output/gff3UcscTest.mapped.gp output/gtfUcscTest.mapped.gp
	${diff} output/gff3UcscTest.unmapped.gp output/gtfUcscTest.unmapped.gp

gff3NcbiTest: mkdirs ${testNcbiLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gff3 data/gencode.v22.annotation.gff3 ${testNcbiLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${gff3ToGenePred} output/$@.unmapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfNcbiTest: mkdirs ${testNcbiLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gtf data/gencode.v22.annotation.gtf ${testNcbiLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/null
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.unmapped.gtf /dev/null
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.unmapped.gtf output/$@.unmapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

# new part naming introduced in V25
gff3ParNamingTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v25.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfParNamingTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v25.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

cmpParNamingTest: gff3ParNamingTest gtfParNamingTest
	${diff} output/gff3ParNamingTest.mapped.gp output/gtfParNamingTest.mapped.gp

# also test without mapinfo or unmapped
gff3UcscSubstituteAuto: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3

# ncRNA and patch subsitution
gff3UcscSubstituteAutoSmallNcRna: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap --useTargetForAutoSmallNonCoding ${targetPatchArg} ${targetGff3Arg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gff3 data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${gff3ToGenePred} output/$@.unmapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gff3UcscSubstituteManOverlap: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap --onlyManualForTargetSubstituteOverlap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGff3Arg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gff3 data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${gff3ToGenePred} output/$@.unmapped.gff3 stdout | ${normalizeGenePred} >output/$@.unmapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfUcscSubstituteManOverlap: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap  --onlyManualForTargetSubstituteOverlap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGtfArg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gtf data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	gtfToGenePred -genePredExt output/$@.unmapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.unmapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.unmapped.gtf output/$@.unmapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

cmpUcscSubstituteManOverlap: gff3UcscSubstituteManOverlap gtfUcscSubstituteManOverlap
	${diff} output/gff3UcscSubstituteManOverlap.mapped.gp output/gtfUcscSubstituteManOverlap.mapped.gp
	${diff} output/gff3UcscSubstituteManOverlap.unmapped.gp output/gtfUcscSubstituteManOverlap.unmapped.gp

# Testing of assigning mapping versions. Use the different results with from NCBI
# to test version numbering.
mappingVerTests: gff3MappingVerBaseTest gtfMappingVerBaseTest \
	gff3MappingVerPrevNoVerTest gtfMappingVerPrevNoVerTest \
	gff3MappingVerPrevDiffMapTest gtfMappingVerPrevDiffMapTest

# initial mapping
gff3MappingVerBaseTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerBaseTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

# tests with previous files don't have mapping version numbers.
# this should produce the same results as gxfMappingVerBaseTest

# remove mapping version number, but don't wack 'appris_alternative_2'
dropMappingVersion = sed -E 's/(=[A-Z]+[0-9]+\.[0-9]+)_[0-9]+;/\1;/g'

gff3MappingVerPrevNoVerTest: mkdirs ${testGencodeLiftOverChains}
	${dropMappingVersion} expected/gff3MappingVerBaseTest.mapped.gff3  >output/$@.prev.gff3
	${gencode_backmap} --oldStyleParIdHack --previousMappedGxf=output/$@.prev.gff3 --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerPrevNoVerTest: mkdirs ${testGencodeLiftOverChains}
	${dropMappingVersion} expected/gtfMappingVerBaseTest.mapped.gtf  >output/$@.prev.gtf
	${gencode_backmap} --oldStyleParIdHack --previousMappedGxf=output/$@.prev.gtf --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info


# Tests with previous files using UCSC chains and new mapping with NCBI. to
# produce different results.  Also add tags to force differences

addTagMatch = /ENSE00001947087/

addTagGff3 = tawk '${addTagMatch}{gsub("tag=","tag=foo,");}{print $$0}'
addTagGtf = tawk '${addTagMatch}{$$0=$$0" tag \"foo\";"}{print $$0}'

gff3MappingVerPrevDiffMapTest: mkdirs ${testNcbiLiftOverChains}
	${addTagGff3} data/gencode.v22.annotation.gff3 > output/$@.input.gff3
	${gencode_backmap} --oldStyleParIdHack --previousMappedGxf=expected/gff3MappingVerBaseTest.mapped.gff3 --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} output/$@.input.gff3 ${testNcbiLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerPrevDiffMapTest: mkdirs ${testNcbiLiftOverChains}
	${addTagGtf} data/gencode.v22.annotation.gtf > output/$@.input.gtf
	${gencode_backmap} --oldStyleParIdHack --previousMappedGxf=expected/gtfMappingVerBaseTest.mapped.gtf --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} output/$@.input.gtf ${testNcbiLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info


# edit to fix names and chrM (don't depend on mkdirs so it doesn't rebuild unless needed)
${testGencodeLiftOverChains}: mkdirs ${testUcscLiftOverChains}
	@mkdir -p output
	../bin/ucscLiftEdit ${testUcscLiftOverChains} data/GCF_000001405.28.assembly.txt data/GCF_000001405.25.assembly.txt ${testGencodeLiftOverChains}

##
# various regression tests
#
##
regressTests: idChangeV25Test gff3V29RegressTest gff3V31RegressTest

##
# Id changes and related in v25
##
# CHD2 - OTT transcript base id moved between genes. This caused chd2 to be dropped due to the target subsitution mechanism though it had already subsituted the gene.
#    OTTHUMT00000415070.1 was in CHD2 in V19/V20 and moved to
#    OTTHUMT00000415070.3 in RP11-437B10.1 in V25
# DUX4 - OTT gene id reused for DUX4L1, Ensemble gene ids changed too.  This causes both
#        failure to map due to non-overlap and failure to substitute because of thinking it's
#        already mapped.
# ZNF577 - manually annotated in V19, no longer manually annotated in V25.
# ENST00000384889.1 move (probably a merge bug)
#   ENST00000384889.1 miRNA in ENSG00000207621.1  V19..V24 MIR224
#   ENST00000384889.1 miRNA in ENSG00000102287.17  V25     GABRE

#        

idChangeV25Test: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap --targetGxf=data/idChange.v19.gff3 ${targetSubstArg} ${headerArg} data/idChange.v25.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info --unmappedGxf=output/$@.unmapped.gff3
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# ENSG00000129932.8/DOHH was not mapped because it it was automatic
# protein-code gene with --useTargetForAutoGenes, ensure it works with
# --useTargetForAutoSmallNonCoding
gff3V29RegressTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap --useTargetForAutoSmallNonCoding ${targetGff3Arg} ${targetSubstArg} ${headerArg} --unmappedGxf=output/$@.unmapped.gff3 data/gencode.v29.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	gff3ToGenePred output/$@.unmapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.unmapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# double mapping of transcript, detected with gff3ToGenePred
gff3V31RegressTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap --useTargetForAutoSmallNonCoding ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v31.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

##
## lift edit
##
ucscLiftEditTest: mkdirs
	../bin/ucscLiftEdit ${testUcscLiftOverChains} data/GCF_000001405.28.assembly.txt data/GCF_000001405.25.assembly.txt output/$@.chain
	${diff} expected/$@.chain output/$@.chain


###
# report program tests
###
reportsTests: mapInfoSumTests

# liftover has more wierd cases, so use for tests
mapInfoSumTests: mapInfoSumGeneAllTest mapInfoSumGeneBiotypeTest mapInfoSumTransBiottypeTest \
	mapInfoSumTransBiocatTest mapInfoSumTransBiocatMultimapTest \
	mapInfoSumTransBiocatTest mapInfoSumTransMatrixGenomicSizeTest \
	mapInfoSumTransBiocatMatrixMapStatusTest \
	mapInfoSumTransBiocatMatRowFreqMapStatusTest \
	mapInfoSumTransBiocatMatColFreqMapStatusTest \
	mapInfoSumSubTargetTest

mapInfoSumGeneAllTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup gene expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumGeneBiotypeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biotype gene expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiottypeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biotype trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiocatTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biocat trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMultimapTest: mkdirs
	../bin/mapInfoSummary --biotypeGroup=biocat --multimapGroup trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransMatrixGenomicSizeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup  --matrixColumnKey=genomicSize trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatrixMapStatusTest: mkdirs
	../bin/mapInfoSummary --matrixColumnKey=mappingStatus --biotypeGroup=biocat trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatRowFreqMapStatusTest: mkdirs
	../bin/mapInfoSummary --rowFrequencies --matrixColumnKey=mappingStatus --biotypeGroup=biocat trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatColFreqMapStatusTest: mkdirs
	../bin/mapInfoSummary --columnFrequencies --matrixColumnKey=mappingStatus --biotypeGroup=biocat trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum

mapInfoSumSubTargetTest: mkdirs
	../bin/mapInfoSummary --targetStatusGroup --substituteTargetGroup trans expected/gff3UcscTest.map-info output/$@.sum
	${diff} expected/$@.sum output/$@.sum




##
# memory leak checks
##
valgrindTest: gff3UcscValgrind gtfUcscValgrind

gff3UcscValgrind: mkdirs
	${valgrind} ${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} data/gencode.v22.annotation.gff3 ${testUcscLiftOverChains} /dev/null /dev/null /dev/null >& output/$@.valgrind.out

gtfUcscValgrind: mkdirs
	${valgrind} ${gencode_backmap} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} data/gencode.v22.annotation.gtf ${testUcscLiftOverChains}  /dev/null /dev/null /dev/null >& output/$@.valgrind.out


liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
v25Dir = /hive/data/genomes/hg38/bed/gencodeV25/data/release_25
v25Gff3 = ${v25Dir}/gencode.v25.chr_patch_hapl_scaff.annotation.gff3.gz
v25Gtf = ${v25Dir}/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz

bigTest: bigGff3UcscTest bigGtfUcscTest

bigGff3UcscTest: mkdirs
	${gencode_backmap} --oldStyleParIdHack --swapMap --unmappedGxf=output/$@.unmapped.gff3 ${v25Gff3} ${liftOverChains} output/$@.mapped.gff3 output/$@.map-info >& output/$@.gff3.log
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${gff3ToGenePred} output/$@.unmapped.gff3 /dev/null

bigGtfUcscTest: mkdirs
	${gencode_backmap} --oldStyleParIdHack --swapMap --unmappedGxf=output/$@.unmapped.gtf ${v25Gtf} ${liftOverChains} output/$@.mapped.gtf output/$@.map-info >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

##
# update test data when annotations change
##
update: mappingChains problemRegions


##
# commands to create new mapping chains when test set is changed
# also creates a PSL from the chain as is useful for debugging
##
mappingChains: ucscMappingChains

ucscMappingChains: tmp/testGrch38.gp
	overlapSelect tmp/testGrch38.gp  /cluster/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz ${testUcscLiftOverChains}
	chainToPslBasic ${testUcscLiftOverChains} tmp/hg38ToHg19.over.psl
	pslSwap tmp/hg38ToHg19.over.psl  tmp/hg19ToHg38.over.psl

# need to only include first ##gff-version header
tmp/testGrch38.gp: ${GRCh38_GFF3s}
	@mkdir -p $(dir $@)
	gff3ToGenePred <(awk 'NR==1||FNR>1{print $$0}' $^) $@

tmp/GRCh37.p13.refseq-gencode.lift: tmp/GCF_000001405.25.assembly.txt
	${ncbiAssemblyReportConvert} --fromIdType=refSeqAccn --toIdType=gencode lift $< $@

tmp/GRCh38.p2.refseq-gencode.lift:  tmp/GCF_000001405.28.assembly.txt
	${ncbiAssemblyReportConvert} --fromIdType=refSeqAccn --toIdType=gencode lift $< $@

tmp/GRCh37.p13.sizes: tmp/GCF_000001405.25.assembly.txt
	${ncbiAssemblyReportConvert} --fromId=refSeqAccn sizes $< $@

tmp/GRCh38.p2.sizes: tmp/GCF_000001405.28.assembly.txt
	${ncbiAssemblyReportConvert} --fromId=refSeqAccn sizes $< $@

tmp/GCF_000001405.25.assembly.txt:
	@mkdir -p tmp
	wget -nv -O $@ ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.assembly.txt

tmp/GCF_000001405.28.assembly.txt:
	wget -nv -O $@ ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.28.assembly.txt


##
# generate problem regions on GRCh37
##
problemRegions: tmp/testGrch37.gp
	wget -nv -O /dev/stdout http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hg19Patch13Patches.txt.gz | zcat | cut -f 2- > tmp/hg19Patch13Patches.bed
	wget -nv -O tmp/Hg19.grcIncidentDb.bb http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/grcIncidentDb/hg19.grcIncidentDb.bb
	bigBedToBed tmp/Hg19.grcIncidentDb.bb stdout | cut -f 1-4 > tmp/Hg19.grcIncidentDb.bed
	(cat tmp/hg19Patch13Patches.bed tmp/Hg19.grcIncidentDb.bed) | overlapSelect -inFmt=bed tmp/testGrch37.gp /dev/stdin data/problemRegions.bed

# need to only include first ##gff-version header
tmp/testGrch37.gp: ${GRCh37_GFF3s}
	@mkdir -p $(dir $@)
	gff3ToGenePred <(awk 'NR==1||FNR>1{print $0}' $^) $@

mkdirs:
	@mkdir -p output

clean:
	rm -rf output tmp
