/*
 * Bring in UCSC browser includes, dealing with C++ conflicts.  Just grab everything needed by any module..
 * Including this first might prevent some conflicts with clang.
 */
#ifndef jkinclude_hh
#define jkinclude_hh
extern "C" {
#define hash jkhash
#include "common.h"
#include "psl.h"
#include "pslTransMap.h"
#include "genomeRangeTree.h"
#include "chain.h"
#include "dnautil.h"
#include "basicBed.h"
#undef hash
}

#endif
