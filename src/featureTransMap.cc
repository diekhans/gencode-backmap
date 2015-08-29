/*
 * Mapping of GxF features using transmap
 */
#include "featureTransMap.hh"
#include "gxf.hh"
#include "pslOps.hh"
#include "pslMapping.hh"

/**
 * Check assumption of feature order being increasing on positive strand and
 * decreasing on negative strand.
 */
bool FeaturesToPsl::checkFeatureOrder(const GxfFeatureVector& features) {
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

/** compute total size of features */
int FeaturesToPsl::sumFeatureSizes(const GxfFeatureVector& features) {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->fEnd - features[i]->fStart)+1;
    }
    return sz;
}

/* build blocks for for alignment. */
void FeaturesToPsl::makePslBlocks(struct psl* psl,
                                  const GxfFeatureVector& features) {
    assert(pslQStrand(psl) == '+');
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < features.size(); iBlk++) {
        const GxfFeature* feature = features[iBlk];
        pslAddBlock(psl, qStart,
                    ((pslTStrand(psl) == '-') ? psl->tStarts[iBlk] = psl->tSize-feature->fEnd : feature->fStart-1),
                    (feature->fEnd - feature->fStart)+1);
        qStart += psl->blockSizes[iBlk];
    }
}

/* construct a PSL */
struct psl* FeaturesToPsl::makeFeaturesPsl(const string& qName,
                                           int qSize, int tStart, int tEnd, int tSize,
                                           const GxfFeatureVector& features) {
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

/* create a psl from a list of features. assumes features are sorter in
 * ascending order.  Return NULL if sequence in features is not in mapping
 * alignments at all. */
struct psl* FeaturesToPsl::toPsl(const string& qName,
                                 int tSize,
                                 const GxfFeatureVector& features) {
    // this does a [1..n] to [0..n) conversion
    if (not checkFeatureOrder(features)) {
        throw invalid_argument("features not in expected order: " + qName);
    }
    int qSize = sumFeatureSizes(features);
    int tStart, tEnd;
    if (features[0]->fStrand == "+") {
        tStart = features[0]->fStart-1;
        tEnd = features[features.size()-1]->fEnd;
    } else {
        tStart = features[features.size()-1]->fStart-1;
        tEnd = features[0]->fEnd;
    }
    
    if (tEnd > tSize) {
        throw invalid_argument("feature " + qName + " end " + toString(tEnd) + " great than " + features[0]->fSeqid
                               + " mapping alignment query size " + toString(tSize) + ", does the mapping alignment need swapped?");
    }
    return makeFeaturesPsl(qName, qSize, tStart, tEnd, tSize, features);
}

/* recursively map a list of PSLs */
void FeatureTransMap::mapPslVector(const PslVector& srcPsls,
                                   int iTransMap,
                                   PslVector& mappedPsls) const {
    for (int i = 0; i < srcPsls.size(); i++) {
        PslVector outPsls = recursiveMapPsl(srcPsls[i], iTransMap);
        mappedPsls.insert(mappedPsls.end(), outPsls.begin(), outPsls.end());
    }
}

/* recursively map a PSL through the coordinate systems */
PslVector FeatureTransMap::recursiveMapPsl(struct psl* srcPsl,
                                           int iTransMap) const {
    PslVector mappedPsls = fTransMaps[iTransMap]->mapPsl(srcPsl);
    if (iTransMap == fTransMaps.size()-1) {
        return mappedPsls;
    } else {
        PslVector mappedPsls2;
        mapPslVector(mappedPsls, iTransMap+1, mappedPsls2);
        mappedPsls.free();
        return mappedPsls2;
    }
    
}

/* Convert features to PSL and map via mapping alignment.  PSL will be kept in
 * the input order of the features, even if this creates a negative target
 * strand. */
PslMapping* FeatureTransMap::mapFeatures(const string& qName,
                                         const GxfFeatureVector& features) const {
    // target is mapping query, which needs to exist to create psl.
    if (not fTransMaps[0]->haveQuerySeq(features[0]->fSeqid)) {
        return NULL;
    }
    int tSize = fTransMaps[0]->getQuerySeqSize(features[0]->fSeqid); // target is mapping query
    struct psl* srcPsl = FeaturesToPsl::toPsl(qName, tSize, features);
    if (srcPsl == NULL) {
         return NULL;
    } else {
        PslVector mappedPsls = recursiveMapPsl(srcPsl, 0);
        return new PslMapping(srcPsl, mappedPsls);
    }
}


/* map a single feature */
PslMapping* FeatureTransMap::mapFeature(const string& qName,
                                        const GxfFeature* feature) const {
    GxfFeatureVector features;
    features.push_back(const_cast<GxfFeature*>(feature));
    return mapFeatures(qName, features);
}
