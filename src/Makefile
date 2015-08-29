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

SRCS = FIOStream.cc gxf.cc gzstream.cc typeOps.cc pslOps.cc frame.cc pslMapping.cc transMap.cc \
	featureTree.cc remapStatus.cc featureTransMap.cc featureMapper.cc geneMapper.cc \
	gencode-backmap.cc

OBJS =  ${SRCS:%.cc=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.cc=%.depend}

all: ${PROG}

${PROG}: ${OBJS}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${KENTLIBS} ${LIBS}

# dependency file is generate as part of compile
${OBJDIR}/%.o: %.cc
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -c -MM -MT $@ $< >$*.depend
	${CXX} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f ${OBJS} ${PROG} ${DEPENDS}
savebak:
	savebak -r ${hgwdev} gencode-backmap Makefile *.cc *.hh HOLD ../tests/data

test:
	cd ../tests && ${MAKE} test

# don't fail on missing dependencies, they are first time the .o is generates
-include ${DEPENDS}
