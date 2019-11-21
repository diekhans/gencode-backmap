#include "remapStatus.hh"
#include <stdexcept>
#include "typeOps.hh"

static const string REMAP_STATUS_NONE_STR = "none";
static const string REMAP_STATUS_FULL_CONTIG_STR = "full_contig";
static const string REMAP_STATUS_FULL_FRAGMENT_STR = "full_fragment";
static const string REMAP_STATUS_PARTIAL_STR = "partial";
static const string REMAP_STATUS_DELETED_STR = "deleted";
static const string REMAP_STATUS_NO_SEQ_MAP_STR = "no_seq_map";
static const string REMAP_STATUS_GENE_CONFLICT_STR = "gene_conflict";
static const string REMAP_STATUS_GENE_SIZE_CHANGE_STR = "gene_size_change";
static const string REMAP_STATUS_AUTO_SMALL_NCRNA_STR = "automatic_small_ncrna_gene";
static const string REMAP_STATUS_AUTOMATIC_GENE_STR = "automatic_gene";
static const string REMAP_STATUS_PSEUDOGENE_STR = "pseudogene";
static const string REMAP_STATUS_INELIGIBLE_STR = "ineligible";
static const string REMAP_STATUS_ERROR_STR = "error";

/* convert a remap status to a string  */
const string& remapStatusToStr(RemapStatus remapStatus) {
    switch (remapStatus) {
        case REMAP_STATUS_NONE:
            return REMAP_STATUS_NONE_STR;
        case REMAP_STATUS_FULL_CONTIG:
            return REMAP_STATUS_FULL_CONTIG_STR;
        case REMAP_STATUS_FULL_FRAGMENT:
            return REMAP_STATUS_FULL_FRAGMENT_STR;
        case REMAP_STATUS_PARTIAL:
            return REMAP_STATUS_PARTIAL_STR;
        case REMAP_STATUS_DELETED:
            return REMAP_STATUS_DELETED_STR;
        case REMAP_STATUS_NO_SEQ_MAP:
            return REMAP_STATUS_NO_SEQ_MAP_STR;
        case REMAP_STATUS_GENE_CONFLICT:
            return REMAP_STATUS_GENE_CONFLICT_STR;
        case REMAP_STATUS_GENE_SIZE_CHANGE:
            return REMAP_STATUS_GENE_SIZE_CHANGE_STR;
        case REMAP_STATUS_AUTO_SMALL_NCRNA:
            return REMAP_STATUS_AUTO_SMALL_NCRNA_STR;
        case REMAP_STATUS_AUTOMATIC_GENE:
            return REMAP_STATUS_AUTOMATIC_GENE_STR;
        case REMAP_STATUS_PSEUDOGENE:
            return REMAP_STATUS_PSEUDOGENE_STR;
        case REMAP_STATUS_INELIGIBLE:
            return REMAP_STATUS_INELIGIBLE_STR;
        case REMAP_STATUS_ERROR:
            return REMAP_STATUS_ERROR_STR;
    }
    return emptyString;
}

/* target status strings */
static const string TARGET_STATUS_NA_STR = "na";
static const string TARGET_STATUS_NEW_STR = "new";
static const string TARGET_STATUS_LOST_STR = "lost";
static const string TARGET_STATUS_OVERLAP_STR = "overlap";
static const string TARGET_STATUS_NONOVERLAP_STR = "nonOverlap";
static const string TARGET_STATUS_ERROR_STR = "error";


/* convert a target status to a string  */
const string& targetStatusToStr(TargetStatus targetStatus) {
    switch (targetStatus) {
        case TARGET_STATUS_NA:
            return TARGET_STATUS_NA_STR;
        case TARGET_STATUS_NEW:
            return TARGET_STATUS_NEW_STR;
        case TARGET_STATUS_LOST:
            return TARGET_STATUS_LOST_STR;
        case TARGET_STATUS_OVERLAP:
            return TARGET_STATUS_OVERLAP_STR;
        case TARGET_STATUS_NONOVERLAP:
            return TARGET_STATUS_NONOVERLAP_STR;
        case TARGET_STATUS_ERROR:
            return TARGET_STATUS_ERROR_STR;
    }
    return emptyString;
}
