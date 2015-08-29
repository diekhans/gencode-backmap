#include "remapStatus.hh"
#include <stdexcept>

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
        case REMAP_STATUS_GENE_EXPAND: {
            static const string status("gene_expand");
            return status;
        }
    }
    static const string notused;
    return notused;
}

