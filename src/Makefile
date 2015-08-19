ROOT = ..

MACH = $(shell uname -m)
SYS = $(shell uname -s)
BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
KENTDIR = ${HOME}/compbio/browser/dev/kent/src
KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACH}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a
LIBS = -lssl -lcrypto -lz
PROG = ${BINDIR}/gencode-backmap

ifeq (${SYS},Darwin)
CXX = g++-mp-4.9
CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare
CXXDEBUG = -g -O0 -fno-default-inline -fno-inline
else
CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
ifeq ($(wildcard ${CXX}),)
CXX = g++
endif
CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare 
CXXDEBUG = -g -gdwarf-3 -O0 -fno-default-inline -fno-inline
LIBS += -lpthread

# autodetect UCSC installation of samtabix:
ifeq (${SAMTABIXDIR},)
    SAMTABIXDIR = /hive/data/outside/samtabix/${MACH}
    ifeq ($(wildcard ${SAMTABIXDIR}),)
        SAMTABIXDIR =
    endif
endif
ifneq (${SAMTABIXDIR},)
    LIBS += ${SAMTABIXDIR}/libsamtabix.a
endif
endif

CXXFLAGS += ${KENTINC} ${CXXDEBUG}

HDRS = FIOStream.hh gxf.hh gzstream.hh typeOps.hh pslOps.hh frame.hh pslMapping.hh transMap.hh \
	gxfFeatureTree.hh remapStatus.hh featureTransMap.hh featureMapper.cc geneMapper.hh \
	jkinclude.hh
SRCS = FIOStream.cc gxf.cc gzstream.cc typeOps.cc pslOps.cc frame.cc pslMapping.cc transMap.cc \
	gxfFeatureTree.cc remapStatus.hh featureTransMap.cc featureMapper.cc geneMapper.cc \
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
	rm -f ${OBJS} ${PROG}
savebak:
	savebak -r ${hgwdev} gencode-backmap Makefile *.cc *.hh HOLD ../tests/data

test:
	cd ../tests && ${MAKE} test
