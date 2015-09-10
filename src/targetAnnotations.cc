#include "targetAnnotations.hh"
#include <stdexcept>

/* load a feature */
void TargetAnnotations::loadFeature(GxfFeature* gxfFeature) {
    string baseId = getBaseId(gxfFeature->getTypeId());
    if (fIdFeatureMap.find(baseId) != fIdFeatureMap.end()) {
        throw invalid_argument("feature with id already exists :" + baseId);
    }
    fIdFeatureMap[baseId] = gxfFeature;
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
        delete it->second;
    }
}

