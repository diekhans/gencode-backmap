ROOT = .

all:
	(cd src && ${MAKE})

test:
	(cd tests && ${MAKE} test)

clean:
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf ${BINDIR}/gencode-backmap.dSYM objs TAGS

etags:
	etags src/*.{hh,cc}
