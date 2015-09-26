#include "targetAnnotations.hh"
#include <stdexcept>

/* load a feature */
void TargetAnnotations::loadFeature(GxfFeature* gxfFeature) {
    string baseId = getBaseId(gxfFeature->getTypeId());
    fIdFeatureMap[baseId].push_back(gxfFeature);
    struct TargetLocationLink* locationLink =  static_cast<struct TargetLocationLink*>(needMem(sizeof(struct TargetLocationLink)));  // zeros memory
    locationLink->feature = gxfFeature;
    genomeRangeTreeAddValList(fLocationMap, toCharStr(gxfFeature->fSeqid),
                              gxfFeature->fStart, gxfFeature->fEnd, locationLink);
}

/* process a record, loading into table or discarding */
void TargetAnnotations::processRecord(GxfRecord* gxfRecord) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        GxfFeature* gxfFeature = dynamic_cast<GxfFeature*>(gxfRecord);
        if ((gxfFeature->fType == GxfFeature::GENE) or (gxfFeature->fType == GxfFeature::TRANSCRIPT)) {
            loadFeature(gxfFeature);
        } else {
            delete gxfRecord;
        }
    } else {
        delete gxfRecord;
    }
}

/* get a target gene or transcript with same base or NULL.
 * special handling for PARs/ */
GxfFeature* TargetAnnotations::get(const string& id,
                                   const string& seqIdForParCheck) const {
    string baseId = getBaseId(id);
    IdFeatureMapConstIter it = fIdFeatureMap.find(baseId);
    if (it == fIdFeatureMap.end()) {
        return NULL;
    } else if (it->second.size() == 2) {
        if (it->second[0]->fSeqid == seqIdForParCheck) {
            return it->second[0];
        } else if (it->second[1]->fSeqid == seqIdForParCheck) {
            return it->second[1];
        } else {
            throw logic_error("PAR target feature hack confused: " + baseId);
        }
    } else {
        return it->second[0];
    }
}

/* find overlapping features */
GxfFeatureVector TargetAnnotations::findOverlappingFeatures(const string& seqid,
                                                            int start,
                                                            int end) {
    GxfFeatureVector overlapping;
    struct range *over, *overs = genomeRangeTreeAllOverlapping(fLocationMap, toCharStr(seqid), start, end);
    for (over = overs; over != NULL; over = over->next) {
        for (struct TargetLocationLink* tll = static_cast<struct TargetLocationLink*>(over->val); tll != NULL; tll = tll->next) {
            overlapping.push_back(tll->feature);
        }
    }
    return overlapping;
}

/* constructor, load gene and transcript objects from a GxF */
TargetAnnotations::TargetAnnotations(const string& gxfFile):
    fLocationMap(genomeRangeTreeNew()) {
    GxfParser gxfParser(gxfFile);
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser.next()) != NULL) {
        processRecord(gxfRecord);
    }
}

/* destructor */
TargetAnnotations::~TargetAnnotations() {
    for (IdFeatureMapIter it = fIdFeatureMap.begin(); it != fIdFeatureMap.end(); it++) {
        it->second.free();
    }

    // free range tree
    struct hashCookie chromCookie = hashFirst(fLocationMap->jkhash);
    struct hashEl *chromEl;
    for (chromEl = hashNext(&chromCookie); chromEl != NULL; chromEl = chromEl->next) {
        struct range *r, *ranges = genomeRangeTreeList(fLocationMap, chromEl->name);
        for (r = ranges; r != NULL; r = r->next) {
            struct TargetLocationLink *tll, *tlls = static_cast<struct TargetLocationLink*>(r->val);
            while ((tll = static_cast<struct TargetLocationLink*>(slPopHead(&tlls))) != NULL) {
                delete tll;
            }
        }
    }
    genomeRangeTreeFree(&fLocationMap);
}

