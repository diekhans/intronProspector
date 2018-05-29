ROOT = ..
include ${ROOT}/defs.mk

all: ${intronProspector}

${intronProspector}: ${OBJS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LIBS}

# dependency file is generate as part of compile
${OBJDIR}/%.o: %.cc
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -c -MM -MT $@ $< >$*.depend
	${CXX} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f ${OBJS} ${PROG} ${DEPENDS}

# Generate manpage.h, which is included to generate help message.  This
# requires pandoc, so it is checked in and only genereated by developers.
doc:  manpage.h


manpage.h: ../docs/intronProspector.1.md
	pandoc -f markdown -t plain  $< >$@.txt.tmp
	awk -f mkManHelp.awk $@.txt.tmp > $@.tmp
	mv $@.tmp $@
	rm -f $@.txt.tmp



# don't fail on missing dependencies, they are first time the .o is generates
-include ${DEPENDS}

# Local Variables:
# mode: makefile-gmake
# End:
