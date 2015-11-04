ROOT = ..
include ${ROOT}/config.mk

targetGxfArg = --targetGxf=data/gencode.v19.annotations.gtf
targetSubstArg = --substituteMissingTargets=V19

valgrind = valgrind --leak-check=full --num-callers=50

testUcscLiftOverChains = data/hg38ToHg19.over.chain
testGencodeLiftOverChains = output/hg38ToHg19.gencode.over.chain

all: test

test: gff3UcscTest gtfUcscTest gff3NcbiTest gtfNcbiTest gff3UcscSkipAutoNcrnaTest\
	mapInfoSumTests ucscLiftEditTest

gff3UcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null
	diff -u expected/$@.mapped.gff3 output/$@.mapped.gff3
	diff -u expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	diff -u expected/$@.map-info output/$@.map-info

gtfUcscTest: mkdirs ${testGencodeLiftOverChains}
	${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gtf ${testGencodeLiftOverChains} output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info
	gtfToGenePred -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/null
	gtfToGenePred -ignoreGroupsWithoutExons output/$@.unmapped.gtf /dev/null
	diff -u expected/$@.mapped.gtf output/$@.mapped.gtf
	diff -u expected/$@.unmapped.gtf output/$@.unmapped.gtf
	diff -u expected/$@.map-info output/$@.map-info

gff3UcscSkipAutoNcrnaTest: mkdirs ${testGencodeLiftOverChains}
	${PROG} --swapMap --skipAutomaticNonCoding ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gff3 ${testGencodeLiftOverChains} output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null
	diff -u expected/$@.mapped.gff3 output/$@.mapped.gff3
	diff -u expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	diff -u expected/$@.map-info output/$@.map-info

# edit to fix names and chrM (don't depend on  mkdirs so it doesn't rebuild unless needed)
${testGencodeLiftOverChains}: ${testUcscLiftOverChains}
	@mkdir -p output
	../bin/ucscLiftEdit ${testUcscLiftOverChains} data/GCF_000001405.28.assembly.txt data/GCF_000001405.25.assembly.txt ${testGencodeLiftOverChains}

gff3NcbiTest: mkdirs
	${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gff3 data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null
	diff -u expected/$@.mapped.gff3 output/$@.mapped.gff3
	diff -u expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	diff -u expected/$@.map-info output/$@.map-info

gtfNcbiTest: mkdirs
	${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gtf data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info
	gtfToGenePred -ignoreGroupsWithoutExons output/$@.mapped.gtf /dev/null
	gtfToGenePred -ignoreGroupsWithoutExons output/$@.unmapped.gtf /dev/null
	diff -u expected/$@.mapped.gtf output/$@.mapped.gtf
	diff -u expected/$@.unmapped.gtf output/$@.unmapped.gtf
	diff -u expected/$@.map-info output/$@.map-info

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
	diff -u expected/$@.sum output/$@.sum

mapInfoSumGeneBiotypeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biotype gene expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiottypeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biotype transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiocatTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup --biotypeGroup=biocat transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMultimapTest: mkdirs
	../bin/mapInfoSummary --biotypeGroup=biocat --multimapGroup transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransMatrixGenomicSizeTest: mkdirs
	../bin/mapInfoSummary --mappingStatusGroup  --matrixColumnKey=genomicSize transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatrixMapStatusTest: mkdirs
	../bin/mapInfoSummary --matrixColumnKey=mappingStatus --biotypeGroup=biocat transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatRowFreqMapStatusTest: mkdirs
	../bin/mapInfoSummary --rowFrequencies --matrixColumnKey=mappingStatus --biotypeGroup=biocat transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumTransBiocatMatColFreqMapStatusTest: mkdirs
	../bin/mapInfoSummary --columnFrequencies --matrixColumnKey=mappingStatus --biotypeGroup=biocat transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

mapInfoSumSubTargetTest: mkdirs
	../bin/mapInfoSummary --targetStatusGroup --substituteTargetGroup transcript expected/gff3UcscTest.map-info output/$@.sum
	diff -u expected/$@.sum output/$@.sum

ucscLiftEditTest: mkdirs
	../bin/ucscLiftEdit ${testUcscLiftOverChains} data/GCF_000001405.28.assembly.txt data/GCF_000001405.25.assembly.txt output/$@.chain
	diff -u expected/$@.chain output/$@.chain

##
# memory leak checks
##
valgrindTest: gff3UcscValgrind gtfUcscValgrind

gff3UcscValgrind: mkdirs
	${valgrind} ${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gff3 ${testUcscLiftOverChains} /dev/null /dev/null /dev/null >& output/$@.valgrind.out

gtfUcscValgrind: mkdirs
	${valgrind} ${gencode_backmap} --swapMap ${targetGxfArg} ${targetSubstArg} data/gencode.v22.annotation.gtf ${testUcscLiftOverChains}  /dev/null /dev/null /dev/null >& output/$@.valgrind.out


liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
v22Dir = /hive/data/genomes/hg38/bed/gencodeV22/data/release_22
v22Gff3 = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gff3.gz
v22Gtf = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz

bigTest: bigGff3UcscTest bigGtfUcscTest bigGff3NcbiTest bigGtfNcbiTest

bigGff3UcscTest: mkdirs
	${gencode_backmap} --swapMap ${v22Gff3} ${liftOverChains} output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info >& output/$@.gff3.log
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

bigGtfUcscTest: mkdirs
	${gencode_backmap} --swapMap ${v22Gtf} ${liftOverChains} output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

bigGff3NcbiTest: mkdirs
	${gencode_backmap} --swapMap ${v22Gff3} data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info >& output/$@.gff3.log
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

bigGtfNcbiTest: mkdirs
	${gencode_backmap} --swapMap ${v22Gtf} data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
