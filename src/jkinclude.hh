/*
 * Bring in UCSC browser includes, dealing with C++ conflicts.  Just grab everything
 */
#ifndef jkinclude_hh
#define jkinclude_hh
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
#include "basicBed.h"
#include "regexHelper.h"
#undef min
#undef max
#undef hash
#undef new
}

#endif
