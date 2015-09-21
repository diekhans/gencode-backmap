#include "remapStatus.hh"
#include <stdexcept>

/* Update a status for a parent from a child.  It starts out as
 * REMAP_STATUS_NONE. */
RemapStatus remapStatusChildUpdate(RemapStatus parentStatus,
                                   RemapStatus childStatus) {

    // Parent status is never changed on these
    switch (parentStatus) {
        case REMAP_STATUS_NO_SEQ_MAP:
        case REMAP_STATUS_GENE_CONFLICT:
        case REMAP_STATUS_GENE_SIZE_CHANGE:
            return parentStatus;
        default:
            break;
    }
    
    // switch statement to generate error if new status is added
    switch (childStatus) {
        case REMAP_STATUS_NONE:
            throw logic_error("child status should not be REMAP_STATUS_NONE");
        case REMAP_STATUS_FULL_CONTIG:
        case REMAP_STATUS_FULL_FRAGMENT:
        case REMAP_STATUS_PARTIAL_CONTIG:
        case REMAP_STATUS_PARTIAL_FRAGMENT:
        case REMAP_STATUS_DELETED:
        case REMAP_STATUS_NO_SEQ_MAP:
        case REMAP_STATUS_GENE_CONFLICT:
        case REMAP_STATUS_GENE_SIZE_CHANGE:
            // worse (highest) wins
            if (parentStatus > childStatus) {
                return parentStatus;
            } else {
                return childStatus;
            }
    }
    throw logic_error("should not make it here");
}

/* convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus) {
    switch (remapStatus) {
        case REMAP_STATUS_NONE: {
            static const string status("none");
            return status;
        }
        case REMAP_STATUS_FULL_CONTIG: {
            static const string status("full_contig");
            return status;
        }
        case REMAP_STATUS_FULL_FRAGMENT: {
            static const string status("full_fragment");
            return status;
        }
        case REMAP_STATUS_PARTIAL_CONTIG: {
            static const string status("partial_contig");
            return status;
        }
        case REMAP_STATUS_PARTIAL_FRAGMENT: {
            static const string status("partial_fragment");
            return status;
        }
        case REMAP_STATUS_DELETED: {
            static const string status("deleted");
            return status;
        }
        case REMAP_STATUS_NO_SEQ_MAP: {
            static const string status("no_seq_map");
            return status;
        }
        case REMAP_STATUS_GENE_CONFLICT: {
            static const string status("gene_conflict");
            return status;
        }
        case REMAP_STATUS_GENE_SIZE_CHANGE: {
            static const string status("gene_size_change");
            return status;
        }
    }
    static const string notused;
    return notused;
}

