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
    REMAP_STATUS_NONE              = 0x0000,  // not set
    REMAP_STATUS_FULL_CONTIG       = 0x0001,  // full remap, contiguous
    REMAP_STATUS_FULL_FRAGMENT     = 0x0002,  // full remap, fragmented
    REMAP_STATUS_PARTIAL           = 0x0004,  // partial remap
    REMAP_STATUS_DELETED           = 0x0008,  // completely deleted
    REMAP_STATUS_NO_SEQ_MAP        = 0x0010,  // no mapping alignments for sequence
    REMAP_STATUS_GENE_CONFLICT     = 0x0020,  // transcripts disagree on sequence/strand
    REMAP_STATUS_GENE_SIZE_CHANGE  = 0x0040,  // change in size of gene has exceeded threshold
    REMAP_STATUS_AUTO_SMALL_NCRNA  = 0x0080,  // automatic small non-coding RNA
    REMAP_STATUS_AUTOMATIC_GENE    = 0x0180,  // automatic gene annotation not mapped
    REMAP_STATUS_PSEUDOGENE        = 0x0200,  // pseudogene annotation not mapped
    REMAP_STATUS_INELIGIBLE        = 0x0400,  // ineligible for mapping
    REMAP_STATUS_NON_PRIMARY       = 0x0900,  // mapped to non-primary chromosome.
    REMAP_STATUS_ERROR             = 0x1000,  // other error
} RemapStatus;

/* Convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus);


/* target status */
typedef enum {
    TARGET_STATUS_NA,         // no target data
    TARGET_STATUS_NEW,        // new gene/transcript
    TARGET_STATUS_LOST,       // gene/transcript not mapped
    TARGET_STATUS_OVERLAP,    // gene/transcript overlaps
    TARGET_STATUS_NONOVERLAP, // gene/transcript doesn't overlaps
    TARGET_STATUS_ERROR       // other error
} TargetStatus;

/* convert a target status to a string  */
const string& targetStatusToStr(TargetStatus targetStatus);

#endif

