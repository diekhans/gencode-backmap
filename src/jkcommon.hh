/*
 * deal with C++ conflicts when include kent common.h
 */
#ifndef jkcommon_hh
#define jkcommon_hh
extern "C" {
#define min jkmin
#define max jkmax
#include "common.h"
#undef min
#undef max
}

#endif
