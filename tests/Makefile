ROOT = ..

BINDIR = ${ROOT}/bin
PROG = ${BINDIR}/gencode-backmap

all: test

test: gff3UcscTest gtfUcscTest

gff3UcscTest: mkdirs
	${PROG} --swapMap data/gencode.v22.annotation.gff3 data/GRCh38-GRCh37.liftover.chains output/$@.mapped.gff3 output/$@.unmapped.gff3
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

gtfUcscTest: mkdirs
	${PROG} --swapMap data/gencode.v22.annotation.gtf data/GRCh38-GRCh37.liftover.chains output/$@.mapped.gtf output/$@.unmapped.gtf
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
v22Dir = /hive/data/genomes/hg38/bed/gencodeV22/data/release_22
v22Gff3 = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gff3.gz
v22Gtf = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz

bigTest: bigGff3UcscTest bigGtfUcscTest

bigGff3UcscTest: mkdirs
	${PROG} --swapMap ${v22Gff3} ${liftOverChains} output/$@.mapped.gff3 output/$@.unmapped.gff3 >& output/$@.gff3.log
	gff3ToGenePred output/$@.mapped.gff3 /dev/null
	gff3ToGenePred output/$@.unmapped.gff3 /dev/null

bigGtfUcscTest: mkdirs
	${PROG} --swapMap ${v22Gtf} ${liftOverChains} output/$@.mapped.gtf output/$@.unmapped.gtf >& output/$@.gtf.log
	gtfToGenePred output/$@.mapped.gtf /dev/null
	gtfToGenePred output/$@.unmapped.gtf /dev/null

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
