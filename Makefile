ROOT = .

pyprogs = mapInfoSummary ncbiAssemblyReportConvert ucscLiftEdit


all:
	(cd src && ${MAKE})

test:
	(cd tests && ${MAKE} test)

clean:
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf ${BINDIR}/gencode-backmap.dSYM objs TAGS

lint:
	${PYTHON} flake8 ${pyprogs:%=bin/%} lib/gencode lib/pycbio


etags:
	etags src/*.{hh,cc}
