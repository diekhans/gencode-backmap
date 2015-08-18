#include "remapStatus.hh"
#include "gxf.hh"

/* Attribute name used for remap status */
const string REMAP_ATTR_NAME = "remap_status";

/* Attribute name used for previous id before remap */
const string REMAP_PREVIOUS_ID = "remap_previous_id";


/* convert a remap status to a AttrVal  */
const AttrVal& remapStatusToAttrVal(RemapStatus remapStatus) {
    switch (remapStatus) {
    case REMAP_STATUS_NONE: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_none");
        return status;
    }
    case REMAP_STATUS_FULL_CONTIG: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_full_contig");
        return status;
    }
    case REMAP_STATUS_FULL_FRAGMENT: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_full_fragment");
        return status;
    }
    case REMAP_STATUS_PARTIAL_CONTIG: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_partial_contig");
        return status;
    }
    case REMAP_STATUS_PARTIAL_FRAGMENT: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_partial_fragment");
        return status;
    }
    case REMAP_STATUS_DELETED: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_deleted");
        return status;
    }
    case REMAP_STATUS_NO_SEQ_MAP: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_no_seq_map");
        return status;
    }
    default:
        throw invalid_argument("BUG: invalid remap status");
    }
}

