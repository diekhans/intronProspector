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

# don't fail on missing dependencies, they are first time the .o is generates
-include ${DEPENDS}

# Local Variables:
# mode: makefile-gmake
# End:
