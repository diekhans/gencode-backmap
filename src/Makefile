ROOT = ..

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
KENTDIR = ${HOME}/compbio/browser/dev/kent/src
KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc

CXXFLAGS += ${KENTINC}

HDRS = FIOStream.hh	gxf.hh		gzstream.hh	typeOps.hh
SRCS = FIOStream.cc	gxf.cc		gzstream.cc	typeOps.cc

OBJS = $(foreach b,$(basename ${SRCS}),${OBJDIR}/${b}.o)

all: ${OBJS}

${OBJDIR}/%.o: %.cc ${HDRS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -c -o $@ $<


savebak:
	savebak -r ${hgwdev} gencode-backmap Makefile *.cc *.hh ../tests/data
