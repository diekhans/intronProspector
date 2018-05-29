ROOT = .
include ${ROOT}/defs.mk

all:
	cd src && ${MAKE} all

test: all
	cd tests && ${MAKE} test

clean:
	cd src && ${MAKE} clean
	cd tests && ${MAKE} clean

savebak:
	savebak -r hgwdev.soe.ucsc.edu intronProspector .

# requires pandoc
doc:
	cd src && ${MAKE} doc
