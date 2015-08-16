ROOT = ..

BINDIR = ${ROOT}/bin
PROG = ${BINDIR}/gencode-backmap

all: test

test: gff3Test gtfTest

gff3Test: mkdirs
	${PROG} data/gencode.v22.annotation.gff3 data/GRCh38-GRCh37.mapping.chains output/$@.mapped.gff3

gtfTest: mkdirs
	${PROG} data/gencode.v22.annotation.gtf data/GRCh38-GRCh37.mapping.chains output/$@.mapped.gtf

liftOverChains = /hive/data/genomes/hg38/bed/liftOver/hg38ToHg19.over.chain.gz
v22Dir = /hive/data/genomes/hg38/bed/gencodeV22/data/release_22
v22Gff3 = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gff3.gz
v22Gtf = ${v22Dir}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz

bigTest: bigGff3Test bigGtfTest

bigGff3Test: mkdirs
	../bin/gencode-backmap --swapMap ${v22Gff3} ${liftOverChains} /dev/null >& output/bigGff3Test.log
bigGtfTest: mkdirs
	../bin/gencode-backmap --swapMap ${v22Gtf} ${liftOverChains} /dev/null >& output/bigGtfTest.log

mkdirs:
	@mkdir -p output

clean:
	rm -f output
