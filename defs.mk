# assumes that ROOT is set to the relative path to directory containing this file.
.PRECIOUS:
OBJDIR = ${ROOT}/objs
BINDIR = ${ROOT}/bin

HTSLIB_PREFIX = /opt/local
HTSLIB_INCL = -I${HTSLIB_PREFIX}/include
HTSLIB_LIB = -L${HTSLIB_PREFIX}/lib -lhts

SRCS = intronProspector.cc junctions_extractor.cc
OBJS =  ${SRCS:%.cc=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.cc=%.depend}
intronProspector = ${BINDIR}/intronProspector

OPT = -g -O0

CXXFLAGS = ${OPT} ${HTSLIB_INCL}
LIBS = ${HTSLIB_LIB}

