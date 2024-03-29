ROOT = ..
include ${ROOT}/config.mk

SHELL = /bin/bash
export SHELLOPTS = pipefail

# to get verbose output: make verbose=--verbose 
verbose = 

### IMPORTANT: test GFF3, if new are added, run make update so chains are rebuilt
GRCh37_GFF3s = data/gencode.v19.annotation.gff3 data/idChange.v19.gff3 \
	data/dup-transcripts.v19.gff3 data/par.V19.gff3 data/transMove.V19.gff3
GRCh38_GFF3s = data/gencode.v22.annotation.gff3 data/idChange.v25.gff3 \
	data/gencode.v29.annotation.gff3 data/gencode.v31.annotation.gff3 \
	data/dup-transcripts.v33.gff3 \
	data/gapped-transcripts.v33.gff3 data/sort.v33.gff3 \
	data/par.V43.gff3 data/par.V44.gff3 data/transMove.V44.gff3

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

gencode_backmap_verb = ${gencode_backmap} ${verbose}

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

test: backmap_tests other_tests
backmap_tests: gff3UcscTest gtfUcscTest cmpUcscTest \
	gff3ParNamingTest gtfParNamingTest cmpParNamingTest \
	gff3NcbiTest gtfNcbiTest \
	gff3UcscSubstituteAuto gff3UcscSubstituteAutoSmallNcRna \
	gff3UcscSubstituteManOverlap gtfUcscSubstituteManOverlap \
	mappingVerTests cmpUcscSubstituteManOverlap \
	gff3V29RegressTest dupTranscriptsV33Test sortStrictV33Test \
	regressTests gapTest parStableIdTests transMoveTests
other_tests: reportsTests ucscLiftEditTest 


gff3UcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info
	${diff} expected/$@.mapped.gp output/$@.mapped.gp

gtfUcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info
	${diff} expected/$@.mapped.gp output/$@.mapped.gp

cmpUcscTest: gff3UcscTest gtfUcscTest
	${diff} output/gff3UcscTest.mapped.gp output/gtfUcscTest.mapped.gp

gff3NcbiTest: mkdirs ${testNcbiLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testNcbiLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfNcbiTest: mkdirs ${testNcbiLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testNcbiLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/null
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

# new part naming introduced in V25
gff3ParNamingTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v25.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfParNamingTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v25.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

cmpParNamingTest: gff3ParNamingTest gtfParNamingTest
	${diff} output/gff3ParNamingTest.mapped.gp output/gtfParNamingTest.mapped.gp

# also test without mapinfo
gff3UcscSubstituteAuto: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3

# ncRNA and patch subsitution
gff3UcscSubstituteAutoSmallNcRna: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap --useTargetForAutoSmallNonCoding ${targetPatchArg} ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gff3UcscSubstituteManOverlap: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap --onlyManualForTargetSubstituteOverlap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info
	${diff} expected/$@.mapped.gp output/$@.mapped.gp

gtfUcscSubstituteManOverlap: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap  --onlyManualForTargetSubstituteOverlap --useTargetForAutoGenes --useTargetForPseudoGenes ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	gtfToGenePred -genePredExt output/$@.mapped.gtf /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info
	${diff} expected/$@.mapped.gp output/$@.mapped.gp

cmpUcscSubstituteManOverlap: gff3UcscSubstituteManOverlap gtfUcscSubstituteManOverlap
	${diff} output/gff3UcscSubstituteManOverlap.mapped.gp output/gtfUcscSubstituteManOverlap.mapped.gp

# Testing of assigning mapping versions. Use the different results with from NCBI
# to test version numbering.
mappingVerTests: gff3MappingVerBaseTest gtfMappingVerBaseTest \
	gff3MappingVerPrevNoVerTest gtfMappingVerPrevNoVerTest \
	gff3MappingVerPrevDiffMapTest gtfMappingVerPrevDiffMapTest

# initial mapping
gff3MappingVerBaseTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerBaseTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --oldStyleParIdHack --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info

# tests with previous files don't have mapping version numbers.
# this should produce the same results as gxfMappingVerBaseTest

# remove mapping version number, but don't wack 'appris_alternative_2'
dropMappingVersion = sed -E 's/(=[A-Z]+[0-9]+\.[0-9]+)_[0-9]+;/\1;/g'

gff3MappingVerPrevNoVerTest: mkdirs ${testGencodeLiftOverChains}
	${dropMappingVersion} expected/gff3MappingVerBaseTest.mapped.gff3  >output/$@.prev.gff3
	${gencode_backmap_verb} --oldStyleParIdHack --previousMappedGxf=output/$@.prev.gff3 --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerPrevNoVerTest: mkdirs ${testGencodeLiftOverChains}
	${dropMappingVersion} expected/gtfMappingVerBaseTest.mapped.gtf  >output/$@.prev.gtf
	${gencode_backmap_verb} --oldStyleParIdHack --previousMappedGxf=output/$@.prev.gtf --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info


# Tests with previous files using UCSC chains and new mapping with NCBI. to
# produce different results.  Also add tags to force differences

addTagMatch = /ENSE00001947087/

addTagGff3 = tawk '${addTagMatch}{gsub("tag=","tag=foo,");}{print $$0}'
addTagGtf = tawk '${addTagMatch}{$$0=$$0" tag \"foo\";"}{print $$0}'

gff3MappingVerPrevDiffMapTest: mkdirs ${testNcbiLiftOverChains}
	${addTagGff3} data/gencode.v22.annotation.gff3 > output/$@.input.gff3
	${gencode_backmap_verb} --oldStyleParIdHack --previousMappedGxf=expected/gff3MappingVerBaseTest.mapped.gff3 --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGff3Arg} ${targetSubstArg} ${headerArg} output/$@.input.gff3 ${testNcbiLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

gtfMappingVerPrevDiffMapTest: mkdirs ${testNcbiLiftOverChains}
	${addTagGtf} data/gencode.v22.annotation.gtf > output/$@.input.gtf
	${gencode_backmap_verb} --oldStyleParIdHack --previousMappedGxf=expected/gtfMappingVerBaseTest.mapped.gtf --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap ${targetGtfArg} ${targetSubstArg} ${headerArg} output/$@.input.gtf ${testNcbiLiftOverChains} output/$@.mapped.gtf output/$@.map-info
	${diff} expected/$@.mapped.gtf output/$@.mapped.gtf
	${diff} expected/$@.map-info output/$@.map-info


# edit to fix names and chrM (don't depend on mkdirs so it doesn't rebuild unless needed)
${testGencodeLiftOverChains}: ${testUcscLiftOverChains}
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
	${gencode_backmap_verb} --swapMap --useTargetForAutoGenes --onlyManualForTargetSubstituteOverlap --targetGxf=data/idChange.v19.gff3 ${targetSubstArg} ${headerArg} data/idChange.v25.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${gff3ToGenePred} output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# ENSG00000129932.8/DOHH was not mapped because it it was automatic
# protein-code gene with --useTargetForAutoGenes, ensure it works with
# --useTargetForAutoSmallNonCoding
gff3V29RegressTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v29.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/stdout | ${normalizeGenePred} > output/$@.mapped.gp
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info
	${diff} expected/$@.mapped.gp output/$@.mapped.gp

# double mapping of transcript, detected with gff3ToGenePred
gff3V31RegressTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding ${targetGff3Arg} ${targetSubstArg} ${headerArg} data/gencode.v31.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# case where transcript moving between genes cause problems
dupTranscriptsV33Test: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding --targetGxf=data/dup-transcripts.v19.gff3 ${targetSubstArg} --previousMappedGxf=data/dup-transcripts.v32lift37.gff3 ${headerArg} data/dup-transcripts.v33.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# std::sort SEGVing because of violating strict ordering
sortStrictV33Test: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap ${headerArg} data/sort.v33.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# gaps in transcripts due to alignment issues
gapTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding ${headerArg} data/gapped-transcripts.v33.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info

# V44 gave unique ids to to chrY PAR, which caused old chrX to be mapped to chrY
parStableIdTests: parStableIdTest

parStableIdTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding --onlyManualForTargetSubstituteOverlap \
            --targetGxf=data/par.V19.gff3 ${targetSubstArg} ${targetPatchArg} --previousMappedGxf=data/par.V43lift37.gff3 ${headerArg} \
	    data/par.V44.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info
	gff3ToGenePred -geneNameAttr=gene_name -rnaNameAttr=transcript_id output/$@.mapped.gff3 output/$@.gp
	genePredCheck -db=hg19 output/$@.gp

transMoveTests: transMoveTest
transMoveTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding --onlyManualForTargetSubstituteOverlap \
            --targetGxf=data/transMove.V19.gff3 ${targetSubstArg} ${targetPatchArg} --previousMappedGxf=data/transMove.V43lift37.gff3 ${headerArg} \
	    data/transMove.V44.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.map-info
	${diff} expected/$@.mapped.gff3 output/$@.mapped.gff3
	${diff} expected/$@.map-info output/$@.map-info
	gff3ToGenePred -geneNameAttr=gene_name -rnaNameAttr=transcript_id output/$@.mapped.gff3 output/$@.gp
	genePredCheck -db=hg19 output/$@.gp



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
	${valgrind} ${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGff3Arg} ${targetSubstArg} data/gencode.v22.annotation.gff3 ${testUcscLiftOverChains} /dev/null /dev/null /dev/null >& output/$@.valgrind.out

gtfUcscValgrind: mkdirs
	${valgrind} ${gencode_backmap_verb} --oldStyleParIdHack --swapMap ${targetGtfArg} ${targetSubstArg} data/gencode.v22.annotation.gtf ${testUcscLiftOverChains}  /dev/null /dev/null /dev/null >& output/$@.valgrind.out


liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
gencodeV19Dir = /hive/data/genomes/hg19/bed/gencodeV19/data/release_19
gencodeV19Gff3 = ${gencodeV19Dir}/gencode.v19.annotation.gff3.gz
bigVer = 44
bigPrevVer = 43
gencodeFullDir = /hive/data/genomes/hg38/bed/gencodeV${bigVer}/hgcImport/data/release_${bigVer}
gencodeFullGff3 = ${gencodeFullDir}/gencode.v${bigVer}.annotation.gff3.gz
gencodePrevLiftDir = /hive/data/genomes/hg19/bed/gencodeV${bigPrevVer}lift37/hgcImport/data/release_${bigPrevVer}lift37
gencodePrevLiftGff3 = ${gencodePrevLiftDir}/gencode.v${bigPrevVer}lift37.annotation.gff3.gz

bigTest: output/update/problemRegions.bed output/update/hg38ToHg19.over.gencode.chain
	@mkdir -p bigout
	${gencode_backmap_verb} --swapMap --useTargetForAutoSmallNonCoding --onlyManualForTargetSubstituteOverlap --substituteMissingTargets=V19 \
	    --targetGxf=${gencodeV19Gff3} --targetPatches=output/update/problemRegions.bed --previousMappedGxf=${gencodePrevLiftGff3} \
	     ${gencodeFullGff3} output/update/hg38ToHg19.over.gencode.chain bigout/gencode.v${bigVer}lift37.annotation.gff3 bigout/gencode.v${bigVer}lift37.map-info.tsv
	${gff3ToGenePred} -geneNameAttr=gene_name -rnaNameAttr=transcript_id bigout/gencode.v${bigVer}lift37.annotation.gff3 bigout/gencode.v${bigVer}lift37.annotation.gp

##
# update test data when annotations change
##
update: mappingChains problemRegions


##
# commands to create new mapping chains when test set is changed
# also creates a PSL from the chain as is useful for debugging
##
mappingChains: output/update/testGrch38.gp
	overlapSelect -inRange output/update/testGrch38.gp /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz ${testUcscLiftOverChains}

##
# generate problem regions on GRCh37
##
problemRegions: output/update/testGrch37.gp output/update/problemRegions.bed
	overlapSelect -inRange output/update/testGrch37.gp output/update/problemRegions.bed data/problemRegions.bed


##
# tmp files used by both update and big targets
##

# Need to convert individually as there are duplicate ids.
# Using /dev/stdout in the loop result in only the last file being 
# the loops stdout, so use UCSC internal mechanism for stdout
output/update/testGrch38.gp: ${GRCh38_GFF3s}
	@mkdir -p $(dir $@)
	for f in $^ ; do gff3ToGenePred $$f stdout ; done > $@.tmp
	mv -f $@.tmp $@

output/update/hg38ToHg19.over.gencode.chain: output/update/GCF_000001405.36_GRCh38.p10_assembly_report.txt output/update/GCF_000001405.13_GRCh37_assembly_report.txt
	@mkdir -p $(dir $@)
	../bin/ucscLiftEdit ${liftOverChains} output/update/GCF_000001405.36_GRCh38.p10_assembly_report.txt output/update/GCF_000001405.13_GRCh37_assembly_report.txt $@

output/update/GCF_000001405.36_GRCh38.p10_assembly_report.txt:
	@mkdir -p $(dir $@)
	wget -nv -o /dev/stderr -O $@.tmp https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/GCF_000001405.36_GRCh38.p10_assembly_report.txt
	mv -f $@.tmp $@

output/update/GCF_000001405.13_GRCh37_assembly_report.txt:
	@mkdir -p $(dir $@)
	wget -nv -o /dev/stderr -O $@.tmp https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt
	mv -f $@.tmp $@

output/update/testGrch37.gp: ${GRCh37_GFF3s}
	@mkdir -p $(dir $@)
	for f in $^ ; do gff3ToGenePred $$f stdout ; done > $@.tmp
	mv -f $@.tmp $@

output/update/problemRegions.bed: output/update/hg19Patch13Patches.bed output/update/Hg19.grcIncidentDb.bed
	cat $^ > $@.tmp
	mv -f $@.tmp $@

output/update/hg19Patch13Patches.bed:
	@mkdir -p $(dir $@)
	wget -nv -o /dev/stderr -O /dev/stdout http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hg19Patch13Patches.txt.gz | zcat | cut -f 2- > $@.tmp
	mv -f $@.tmp $@

output/update/Hg19.grcIncidentDb.bed:
	@mkdir -p $(dir $@)
	wget -nv -o /dev/stderr -O output/update/Hg19.grcIncidentDb.bb.tmp http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/grcIncidentDb/hg19.grcIncidentDb.bb
	mv -f output/update/Hg19.grcIncidentDb.bb.tmp output/update/Hg19.grcIncidentDb.bb
	bigBedToBed output/update/Hg19.grcIncidentDb.bb stdout | cut -f 1-4 > $@.tmp
	mv -f $@.tmp $@


mkdirs:
	@mkdir -p output

clean:
	rm -rf output tmp bigout
