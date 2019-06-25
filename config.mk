MACH = $(shell uname -m)
SYS = $(shell uname -s)
SHELL = /bin/bash

PYTHON = python3

# edit to set to UCSC browser kent/src
KENTDIR = ${HOME}/kent/src

KENTINC = -I${KENTDIR}/inc -I${KENTDIR}/hg/inc
KENTLIBDIR = ${KENTDIR}/lib/${MACH}
KENTLIBS = ${KENTLIBDIR}/jkhgap.a ${KENTLIBDIR}/jkweb.a ${KENTDIR}/htslib/libhts.a
LIBS = -lssl -lcrypto -lz -lpthread

ifeq (${SYS},Darwin)
    USE_CLANG = yes
else
    USE_CLANG = no
endif

ifeq (${USE_CLANG}, yes)
    CXX = clang++
    CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare -Wno-format-security
    CXXDEBUG = -g -O0 -fno-inline
else
    # gcc
    ifeq (${SYS},Darwin)
        CXX = g++-mp-8
    else ifneq ($(wildcard /opt/rh/devtoolset-2/root/usr/bin/g++),)
        CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
    else
        CXX = g++
    endif
    CXXFLAGS = -std=c++11 -Wall -Werror -Wno-sign-compare 
    CXXDEBUG = -g -gdwarf-3 -O0 -fno-default-inline -fno-inline
endif
#CXXDEBUG += -pg

CXXFLAGS += ${KENTINC} ${CXXDEBUG}

BINDIR = ${ROOT}/bin
OBJDIR = ${ROOT}/objs
gencode_backmap = ${BINDIR}/gencode-backmap
gencodeAttrsStats = ${BINDIR}/gencodeAttrsStats
