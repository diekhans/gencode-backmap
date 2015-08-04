/*
 * Bring in UCSC browser includes, dealing with C++ conflicts.  Just grab everything
 */
#ifndef jkcommon_hh
#define jkcommon_hh
extern "C" {
#define min jkmin
#define max jkmax
#define hash jkhash
#define new jknew
#include "common.h"
#include "psl.h"
#include "pslTransMap.h"
#include "genomeRangeTree.h"
#include "chain.h"
#include "dnautil.h"
#undef min
#undef max
#undef hash
#undef new
}

#endif
