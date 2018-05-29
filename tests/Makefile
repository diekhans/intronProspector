ROOT = ..
include ${ROOT}/defs.mk

diff = diff -u

all:

test: test1

test1: mkout
	${intronProspector} --junction-bed=output/$@.juncs.bed --intron-bed=output/$@.introns.bed input/test_basic.sam 
	${diff} expected/$@.juncs.bed output/$@.juncs.bed
	${diff} expected/$@.introns.bed output/$@.introns.bed

mkout:
	@mkdir -p output

clean:
	rm -rf output

