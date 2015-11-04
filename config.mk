MACH = $(shell uname -m)
SYS = $(shell uname -s)

# edit to set to UCSC browser kent/src
KENTDIR = ${HOME}/kent/src

KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACH}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a
LIBS = -lssl -lcrypto -lz -lpthread

# autodetect UCSC installation of samtabix
ifeq (${SAMTABIXDIR},)
    SAMTABIXDIR = /hive/data/outside/samtabix/${MACH}
    ifeq ($(wildcard ${SAMTABIXDIR}),)
        SAMTABIXDIR =
    endif
endif
ifneq (${SAMTABIXDIR},)
    LIBS += ${SAMTABIXDIR}/libsamtabix.a
endif


ifeq (${SYS},Darwin)
    CXX = g++-mp-4.9
    CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare
    CXXDEBUG = -g -O0 -fno-default-inline -fno-inline
else ifneq ($(wildcard /opt/rh/devtoolset-2/root/usr/bin/g++),)
    CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
else
    CXX = g++
endif
CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare 
CXXDEBUG = -g -gdwarf-3 -O0 -fno-default-inline -fno-inline

CXXFLAGS += ${KENTINC} ${CXXDEBUG}

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
gencode_backmap = ${BINDIR}/gencode-backmap
