ROOT = .

pyprogs = mapInfoSummary ncbiAssemblyReportConvert ucscLiftEdit


all:
	(cd src && ${MAKE})

test:
	(cd tests && ${MAKE} test)

clean:
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf bin/gencode-backmap.dSYM objs TAGS
	rm -rf $$(find lib -name __pycache__)
lint:
	${PYTHON} flake8 ${pyprogs:%=bin/%} lib/gencode lib/pycbio_local


etags:
	etags src/*.{hh,cc}
