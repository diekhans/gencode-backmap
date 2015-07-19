ROOT = ..

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
KENTDIR = ${HOME}/compbio/browser/dev/kent/src
KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
MACH = $(shell uname -m)
KENTLIBDIR = ${KENTDIR}/lib/${MACH}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a
LIBS = -lssl -lcrypto -lz
PROG = ${BINDIR}/gencode-backmap

CXX = g++-mp-4.8 -std=c++11
CXXDEBUG = -g -O0 -fno-default-inline -fno-inline

CXX = c++
CXXDEBUG = -g -O0 -fno-inline

CXXFLAGS += ${KENTINC} -Wall -Werror ${CXXDEBUG}

HDRS = FIOStream.hh gxf.hh  gzstream.hh typeOps.hh transMapper.hh geneMapper.hh
SRCS = FIOStream.cc gxf.cc  gzstream.cc typeOps.cc transMapper.cc geneMapper.cc \
	gencode-backmap.cc

OBJS = $(foreach b,$(basename ${SRCS}),${OBJDIR}/${b}.o)

all: ${PROG}

${PROG}: ${OBJS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${KENTLIBS} ${LIBS}

${OBJDIR}/%.o: %.cc ${HDRS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f ${OBJS}

savebak:
	savebak -r ${hgwdev} gencode-backmap Makefile *.cc *.hh ../tests/data
