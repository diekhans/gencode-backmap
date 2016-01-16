#include "bedMap.hh"
#include "jkinclude.hh"
#include "typeOps.hh"

/* constructor */
BedMap::BedMap(const string& bedFile):
    fLocationMap(genomeRangeTreeNew()) {
    struct bed* allRecs = bedLoadNAll(toCharStr(bedFile), 3);
    struct bed* rec;
    while ((rec = static_cast<struct bed*>(slPopHead(&allRecs))) != NULL) {
        genomeRangeTreeAddValList(fLocationMap, rec->chrom, rec->chromStart, rec->chromEnd, rec);
    }
}

/* destructor */
BedMap::~BedMap() {
    struct hashCookie chromCookie = hashFirst(fLocationMap->jkhash);
    struct hashEl *chromEl;
    while ((chromEl = hashNext(&chromCookie)) != NULL) {
        struct range *ranges = genomeRangeTreeList(fLocationMap, chromEl->name);
        for (struct range *r = ranges; r != NULL; r = r->next) {
            struct bed* beds = static_cast<struct bed*>(r->val);
            bedFreeList(&beds);
        }
    }
    genomeRangeTreeFree(&fLocationMap);
}

/* check for overlap on a range */
bool BedMap::anyOverlap(const string& seqid,
                        int start,
                        int end) const {
    struct range *overs = genomeRangeTreeAllOverlapping(fLocationMap, toCharStr(seqid), start, end);
    return (overs != NULL);
}
