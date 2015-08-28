#include "remapStatus.hh"
#include <stdexcept>

/* Attribute name used for remap status */
const string REMAP_ATTR_NAME = "remap_status";

/* Attribute name used for original id before remap */
const string REMAP_ORIGINAL_ID = "remap_original_id";

/* convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus) {
    switch (remapStatus) {
    case REMAP_STATUS_NONE: {
        static const string status("remap_status_none");
        return status;
    }
    case REMAP_STATUS_FULL_CONTIG: {
        static const string status("remap_status_full_contig");
        return status;
    }
    case REMAP_STATUS_FULL_FRAGMENT: {
        static const string status("remap_status_full_fragment");
        return status;
    }
    case REMAP_STATUS_PARTIAL_CONTIG: {
        static const string status("remap_status_partial_contig");
        return status;
    }
    case REMAP_STATUS_PARTIAL_FRAGMENT: {
        static const string status("remap_status_partial_fragment");
        return status;
    }
    case REMAP_STATUS_DELETED: {
        static const string status("remap_status_deleted");
        return status;
    }
    case REMAP_STATUS_NO_SEQ_MAP: {
        static const string status("remap_status_no_seq_map");
        return status;
    }
    default:
        throw invalid_argument("BUG: invalid remap status");
    }
}

