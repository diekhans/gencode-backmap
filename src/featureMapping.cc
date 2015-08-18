/*
 * Mapping of features
 */
#include "featureMapping.hh"
#include <iostream>
#include <algorithm> 

/* compute the remap status of the feature */
RemapStatus FeatureMapping::calcRemapStatus() const {
    if (not fSrcSeqInMapping) {
        // couldn't even try mapping, chrom not in map
        return REMAP_STATUS_NO_SEQ_MAP;
    } else if (fMappedFeatures.size() == 0) {
        assert(fUnmappedFeatures.size() > 0);
        // nothing mapped
        return REMAP_STATUS_DELETED;
    } else if (fUnmappedFeatures.size() == 0) {
        // full mapped
        if (fMappedFeatures.size() == 1) {
            return REMAP_STATUS_FULL_CONTIG;
        } else {
            return REMAP_STATUS_FULL_FRAGMENT;
        }
    } else {
        // partially mapped
        if (fMappedFeatures.size() == 1) {
            return REMAP_STATUS_PARTIAL_CONTIG;
        } else {
            return REMAP_STATUS_PARTIAL_FRAGMENT;
        }
    }
}

/* print for debugging */
void FeatureMapping::dump(ostream& fh) const {
    const string status = remapStatusToStr(calcRemapStatus());
    // combine and sort into order
    GxfFeatureVector both;
    both.insert(both.end(), fMappedFeatures.begin(), fMappedFeatures.end());
    both.insert(both.end(), fUnmappedFeatures.begin(), fUnmappedFeatures.end());
    for (int i = 0; i < both.size(); i++) {
        bool isMapped = std::find(fMappedFeatures.begin(), fMappedFeatures.end(), both[i]) != fMappedFeatures.end();
        fh << (isMapped ? "mapped" : "unmapped") << "\t" << status
               << "\t" << both[i]->toString() << endl;
    }
}

/* print for debugging */
void FeatureMappingSet::dump(ostream& fh) const {
    for (int i = 0; i < size(); i++) {
        (*this)[i]->dump(fh);
    }
}
