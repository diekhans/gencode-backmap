/*
 * Mapping of GxF features using transmap
 */
#include "featureTransMap.hh"
#include "gxf.hh"
#include "pslOps.hh"

/** compute total size of features */
int FeatureTransMap::sumFeatureSizes(const GxfFeatureVector& features) const {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->fEnd - features[i]->fStart)+1;
    }
    return sz;
}

/* build blocks for for alignment. */
void FeatureTransMap::makePslBlocks(struct psl* psl,
                                    const GxfFeatureVector& features) const {
    assert(pslQStrand(psl) == '+');
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < features.size(); iBlk++) {
        const GxfFeature* feature = features[iBlk];
        psl->blockSizes[iBlk] = (feature->fEnd - feature->fStart)+1;
        psl->qStarts[iBlk] = qStart;
        if (pslTStrand(psl) == '-') {
            psl->tStarts[iBlk] = psl->tSize-feature->fEnd;
        } else {
            psl->tStarts[iBlk] = (feature->fStart-1);
        }
        psl->match += psl->blockSizes[iBlk];
        qStart += psl->blockSizes[iBlk];
        psl->blockCount++;
    }
}


/* create a psl from a list of features. assumes features are sorter in
 * ascending order */
struct psl* FeatureTransMap::featuresToPsl(const string& qName,
                                           const GxfFeatureVector& features) const {
    // this does a [1..n] to [0..n) conversion
    int qSize = sumFeatureSizes(features);
    int tStart, tEnd;
    if (features[0]->fStrand == "+") {
        tStart = features[0]->fStart-1;
        tEnd = features[features.size()-1]->fEnd;
    } else {
        tStart = features[features.size()-1]->fStart-1;
        tEnd = features[0]->fEnd;
    }
    
    int tSize = fTransMap->getQuerySeqSize(features[0]->fSeqid); // target is mapping query
    char strand[3] = {'+', features[0]->fStrand[0], '\0'};

    struct psl* psl = pslNew(toCharStr(qName), qSize, 0, qSize,
                             toCharStr(features[0]->fSeqid), tSize, tStart, tEnd,
                             strand, features.size(), 0);
    makePslBlocks(psl, features);
    if (pslCheck(toCharStr("converted GxF"), stderr, psl) > 0) {
        throw invalid_argument("invalid PSL created: " + pslToString(psl));
    }
    return psl;
}

/**
 * Check assumption of feature order being increasing on positive strand and
 * decreasing on negative strand.
 */
bool FeatureTransMap::checkFeatureOrder(const GxfFeatureVector& features) const {
    for (int iFeat = 1; iFeat < features.size(); iFeat++) {
        if (features[0]->fStrand == "+") {
            if (features[iFeat]->fStart <= features[iFeat-1]->fEnd) {
                return false;
            }
        } else {
            if (features[iFeat]->fEnd >= features[iFeat-1]->fStart) {
                return false;
            }
        }
    }
    return true;
}

/* Map features to list of PSLs.  PSL will be kept in the input order
 * of the features, even if this creates a negative target strand. */
PslMapping* FeatureTransMap::mapFeatures(const string& qName,
                                         const GxfFeatureVector& features) const {
    if (not checkFeatureOrder(features)) {
        throw invalid_argument("features not in expected order: " + qName);
    }
    
    struct psl* inPsl = featuresToPsl(qName, features);
    return fTransMap->mapPsl(inPsl);
}
