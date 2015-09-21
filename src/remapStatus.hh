/*
 * status of remapping.
 */
#ifndef remapStatus_hh
#define remapStatus_hh
#include <string>
using namespace std;

/* status of remap of feature, these are in order from best to worst after
 * NONE.  This is used to decide status from children.
 */
typedef enum {
    REMAP_STATUS_NONE,               // not set
    REMAP_STATUS_FULL_CONTIG,        // full remap, contiguous
    REMAP_STATUS_FULL_FRAGMENT,      // full remap, fragmented
    REMAP_STATUS_PARTIAL_CONTIG,     // partial remap, contiguous
    REMAP_STATUS_PARTIAL_FRAGMENT,   // partial remap, fragmented
    REMAP_STATUS_DELETED,            // completely deleted
    REMAP_STATUS_NO_SEQ_MAP,         // no mapping alignments for sequence
    REMAP_STATUS_GENE_CONFLICT,      // transcripts disagree on sequence/strand
    REMAP_STATUS_GENE_SIZE_CHANGE,    // change in size of gene has exceeded threshold
} RemapStatus;

/* Update a status for a parent from a child.  It starts out as
 * REMAP_STATUS_NONE. */
RemapStatus remapStatusChildUpdate(RemapStatus parentStatus,
                                   RemapStatus childStatus);

/* convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus);
   

#endif

