/*
 * status of remapping.
 */
#ifndef remapStatus_hh
#define remapStatus_hh
#include <string>
using namespace std;

/* status of remap of feature */
typedef enum {
    REMAP_STATUS_NONE,               // not set
    REMAP_STATUS_FULL_CONTIG,        // full remap, contiguous
    REMAP_STATUS_FULL_FRAGMENT,      // full remap, fragmented
    REMAP_STATUS_PARTIAL_CONTIG,     // partial remap, contiguous
    REMAP_STATUS_PARTIAL_FRAGMENT,   // partial remap, fragmented
    REMAP_STATUS_DELETED,            // completely deleted
    REMAP_STATUS_NO_SEQ_MAP          // no mapping alignments for sequence
} RemapStatus;

/* Attribute name used for remap status */
extern const string REMAP_ATTR_NAME;

/* Attribute name used for previous id before remap */
extern const string REMAP_PREVIOUS_ID;

/* convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus);
   

#endif

