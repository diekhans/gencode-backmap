/*
 * status of remapping.
 */
#ifndef remapStatus_hh
#define remapStatus_hh
#include <string>
using namespace std;

/* Status of remap of feature.
 */
typedef enum {
    REMAP_STATUS_NONE              = 0x00,  // not set
    REMAP_STATUS_FULL_CONTIG       = 0x01,  // full remap, contiguous
    REMAP_STATUS_FULL_FRAGMENT     = 0x02,  // full remap, fragmented
    REMAP_STATUS_PARTIAL           = 0x04,  // partial remap
    REMAP_STATUS_DELETED           = 0x08,  // completely deleted
    REMAP_STATUS_NO_SEQ_MAP        = 0x10,  // no mapping alignments for sequence
    REMAP_STATUS_GENE_CONFLICT     = 0x20,  // transcripts disagree on sequence/strand
    REMAP_STATUS_GENE_SIZE_CHANGE  = 0x40,  // change in size of gene has exceeded threshold
    REMAP_STATUS_AUTOMATIC_GENE    = 0x80   // small, automatic non-coding RNA
} RemapStatus;

/* Convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus);


/* target status */
typedef enum {
    TARGET_STATUS_NA,         // no target data
    TARGET_STATUS_NEW,        // new gene/transcript
    TARGET_STATUS_LOST,       // gene/transcript not mapped
    TARGET_STATUS_OVERLAP,    // gene/transcript overlaps
    TARGET_STATUS_NONOVERLAP  // gene/transcript doesn't overlaps
} TargetStatus;

/* convert a target status to a string  */
const string& targetStatusToStr(TargetStatus targetStatus);

#endif

