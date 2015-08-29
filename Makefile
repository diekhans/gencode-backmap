ROOT = .

all:
	(cd src && ${MAKE})

test:
	(cd tests && ${MAKE} test)

clean:
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
