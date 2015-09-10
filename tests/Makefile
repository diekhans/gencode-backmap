ROOT = ..

BINDIR = ${ROOT}/bin
PROG = ${BINDIR}/gencode-backmap

targetGxfArg = --targetGxf=data/gencode.v19.annotations.gtf

valgrind = valgrind --leak-check=full --num-callers=50

all: test

test: gff3UcscTest gtfUcscTest gff3NcbiTest gtfNcbiTest mapInfoSumTest

gff3UcscTest: mkdirs
	${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gff3 data/GRCh38-GRCh37.liftover.chain output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null
	diff -u expected/$@.mapped.gff3 output/$@.mapped.gff3
	diff -u expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	diff -u expected/$@.map-info output/$@.map-info

gtfUcscTest: mkdirs
	${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gtf data/GRCh38-GRCh37.liftover.chain output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null
	diff -u expected/$@.mapped.gtf output/$@.mapped.gtf
	diff -u expected/$@.unmapped.gtf output/$@.unmapped.gtf
	diff -u expected/$@.map-info output/$@.map-info

gff3NcbiTest: mkdirs
	${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gff3 data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null
	diff -u expected/$@.mapped.gff3 output/$@.mapped.gff3
	diff -u expected/$@.unmapped.gff3 output/$@.unmapped.gff3
	diff -u expected/$@.map-info output/$@.map-info

gtfNcbiTest: mkdirs
	${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gtf data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null
	diff -u expected/$@.mapped.gtf output/$@.mapped.gtf
	diff -u expected/$@.unmapped.gtf output/$@.unmapped.gtf
	diff -u expected/$@.map-info output/$@.map-info

# liftover has more wierd cases, so use for tests
mapInfoSumTest: mkdirs
	../bin/mapInfoSummary expected/gff3UcscTest.map-info output/$@.map-info-summary
	diff -u expected/$@.map-info-summary output/$@.map-info-summary

valgrindTest: gff3UcscValgrind gtfUcscValgrind

gff3UcscValgrind: mkdirs
	${valgrind} ${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gff3 data/GRCh38-GRCh37.liftover.chain /dev/null /dev/null /dev/null >& output/$@.valgrind.out

gtfUcscValgrind: mkdirs
	${valgrind} ${PROG} --swapMap ${targetGxfArg} data/gencode.v22.annotation.gtf data/GRCh38-GRCh37.liftover.chain  /dev/null /dev/null /dev/null >& output/$@.valgrind.out


liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
v22Dir = /hive/data/genomes/hg38/bed/gencodeV22/data/release_22
v22Gff3 = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gff3.gz
v22Gtf = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz

bigTest: bigGff3UcscTest bigGtfUcscTest bigGff3NcbiTest bigGtfNcbiTest

bigGff3UcscTest: mkdirs
	${PROG} --swapMap ${v22Gff3} ${liftOverChains} output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info >& output/$@.gff3.log
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

bigGtfUcscTest: mkdirs
	${PROG} --swapMap ${v22Gtf} ${liftOverChains} output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

bigGff3NcbiTest: mkdirs
	${PROG} --swapMap ${v22Gff3} data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gff3 output/$@.unmapped.gff3 output/$@.map-info >& output/$@.gff3.log
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

bigGtfNcbiTest: mkdirs
	${PROG} --swapMap ${v22Gtf} data/GRCh38.p2-GRCh37.p13.gencode.psl output/$@.mapped.gtf output/$@.unmapped.gtf output/$@.map-info >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
