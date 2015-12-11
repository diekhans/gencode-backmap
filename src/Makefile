ROOT = ..
include ${ROOT}/config.mk



SRCS = FIOStream.cc gxf.cc gzstream.cc typeOps.cc pslOps.cc frame.cc pslMapping.cc transMap.cc \
	featureTree.cc remapStatus.cc  annotationSet.cc featureTransMap.cc \
	featureMapper.cc transcriptMapper.cc geneMapper.cc bedMap.cc gencode-backmap.cc

OBJS =  ${SRCS:%.cc=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.cc=%.depend}

all: ${gencode_backmap}

${gencode_backmap}: ${OBJS}
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
	savebak -r ${hgwdev} gencode-backmap Makefile *.cc *.hh ../tests/data

test:
	cd ../tests && ${MAKE} test

# don't fail on missing dependencies, they are first time the .o is generates
-include ${DEPENDS}


# Local Variables:
# mode: makefile-gmake
# End:
