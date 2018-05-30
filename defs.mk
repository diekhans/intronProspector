# assumes that ROOT is set to the relative path to directory containing this file.
.PRECIOUS:
OBJDIR = ${ROOT}/objs
BINDIR = ${ROOT}/bin

ifneq (${MED_OPT},)
   HTSLIB_PREFIX = ${MED_OPT}
   CXX =  /opt/rh/devtoolset-6/root/usr/bin/g++ 
else
   HTSLIB_PREFIX = /opt/local
   HTSLIB_PREFIX = /Users/markd/compbio/gencode/projs/icedb/array-express/src/local
endif
HTSLIB_INCL = -I${HTSLIB_PREFIX}/include
HTSLIB_LIB = -L${HTSLIB_PREFIX}/lib -lhts

SRCS = intronProspector.cc junctions_extractor.cc
OBJS =  ${SRCS:%.cc=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.cc=%.depend}
intronProspector = ${BINDIR}/intronProspector

OPT = -g -O0 -Wall

CXXFLAGS = ${OPT} ${HTSLIB_INCL}
LIBS = ${HTSLIB_LIB}

