ROOT = ..
include ${ROOT}/defs.mk

# Test run command in various ways and compares to expected output

diff = diff -u

all:

basic_sam = input/test_basic.sam
basic_bam = output/test_basic.bam

test: test1 test2 test3

# basic unstranded
test1: mkout
	${intronProspector} --junction-bed=output/$@.juncs.bed --intron-bed=output/$@.introns.bed --pass-through=output/$@.sam ${basic_sam}
	${diff} expected/$@.juncs.bed output/$@.juncs.bed
	${diff} expected/$@.introns.bed output/$@.introns.bed
	${diff} output/test1.sam output/$@.sam

# RF strandness
test2: mkout
	${intronProspector} --strandness=RF --junction-bed=output/$@.juncs.bed --intron-bed=output/$@.introns.bed ${basic_sam}
	${diff} expected/$@.juncs.bed output/$@.juncs.bed
	${diff} expected/$@.introns.bed output/$@.introns.bed

# test1 with BAM
test3: mkout ${basic_bam}
	${intronProspector} --junction-bed=output/$@.juncs.bed --intron-bed=output/$@.introns.bed --pass-through=output/$@.bam ${basic_bam}
	${diff} expected/test1.juncs.bed output/$@.juncs.bed
	${diff} expected/test1.introns.bed output/$@.introns.bed
	samtools view -h output/$@.bam >output/$@.sam
	${diff} output/test1.sam output/$@.sam

mkout:
	@mkdir -p output

output/%.bam: input/%.sam
	@mkdir -p output
	samtools view -b $< >$@.tmp
	mv -f $@.tmp $@

clean:
	rm -rf output

