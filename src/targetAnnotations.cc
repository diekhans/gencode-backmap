#include "targetAnnotations.hh"
#include <stdexcept>

/* load a feature */
void TargetAnnotations::loadFeature(GxfFeature* gxfFeature) {
    string baseId = getBaseId(gxfFeature->getTypeId());
    fIdFeatureMap[baseId].push_back(gxfFeature);
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

/* constructor, load gene and transcript objects from a GxF */
TargetAnnotations::TargetAnnotations(const string& gxfFile) {
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

}

