ROOT = ..

BINDIR = ${ROOT}/bin
PROG = ${BINDIR}/gencode-backmap

all: test

test: gff3Test gtfTest

gff3Test: mkdirs
	${PROG} data/gencode.v22.annotation.gff3 data/GRCh38-GRCh37.mapping.chains output/$@.mapped.gff3

gtfTest: mkdirs
	${PROG} data/gencode.v22.annotation.gtf data/GRCh38-GRCh37.mapping.chains output/$@.mapped.gtf

mkdirs:
	@mkdir -p output

clean:
	rm -f output
