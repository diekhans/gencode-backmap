/*
 * Mapping of GxF features using transmap
 */
#include "featureTransMap.hh"
#include "feature.hh"
#include "pslOps.hh"
#include "pslMapping.hh"

/**
 * Check assumption of feature order being increasing on positive strand and
 * decreasing on negative strand.
 */
bool FeaturesToPsl::checkFeatureOrder(const FeatureVector& features) {
    for (int iFeat = 1; iFeat < features.size(); iFeat++) {
        if (features[0]->getStrand() == "+") {
            if (features[iFeat]->getStart() <= features[iFeat-1]->getEnd()) {
                return false;
            }
        } else {
            if (features[iFeat]->getEnd() >= features[iFeat-1]->getStart()) {
                return false;
            }
        }
    }
    return true;
}

/** compute total size of features */
int FeaturesToPsl::sumFeatureSizes(const FeatureVector& features) {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->getEnd() - features[i]->getStart())+1;
    }
    return sz;
}

/* build blocks for for alignment. */
void FeaturesToPsl::makePslBlocks(struct psl* psl,
                                  const FeatureVector& features) {
    assert(pslQStrand(psl) == '+');
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < features.size(); iBlk++) {
        const Feature* feature = features[iBlk];
        pslAddBlock(psl, qStart,
                    ((pslTStrand(psl) == '-') ? psl->tStarts[iBlk] = psl->tSize-feature->getEnd() : feature->getStart()-1),
                    (feature->getEnd() - feature->getStart())+1);
        qStart += psl->blockSizes[iBlk];
    }
}

/* construct a PSL */
struct psl* FeaturesToPsl::makeFeaturesPsl(const string& qName,
                                           int qSize, int tStart, int tEnd, int tSize,
                                           const FeatureVector& features) {
    char strand[3] = {'+', features[0]->getStrand()[0], '\0'};

    struct psl* psl = pslNew(toCharStr(qName), qSize, 0, qSize,
                             toCharStr(features[0]->getSeqid()), tSize, tStart, tEnd,
                             strand, features.size(), 0);
    makePslBlocks(psl, features);
    if (pslCheck(toCharStr("converted GxF"), stderr, psl) > 0) {
        throw invalid_argument("invalid PSL created: " + pslToString(psl));
    }
    return psl;
}

/* create a psl from a list of features. Assumes features are sorter in
 * ascending order.  */
struct psl* FeaturesToPsl::toPsl(const string& qName,
                                 int tSize,
                                 const FeatureVector& features) {
    // this does a [1..n] to [0..n) conversion
    if (not checkFeatureOrder(features)) {
        throw invalid_argument("features not in expected order: " + qName);
    }
    int qSize = sumFeatureSizes(features);
    int tStart, tEnd;
    if (features[0]->getStrand() == "+") {
        tStart = features[0]->getStart()-1;
        tEnd = features[features.size()-1]->getEnd();
    } else {
        tStart = features[features.size()-1]->getStart()-1;
        tEnd = features[0]->getEnd();
    }
    
    if (tEnd > tSize) {
        throw invalid_argument("feature " + qName + " end " + toString(tEnd) + " great than " + features[0]->getSeqid()
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
                                         const FeatureVector& features) const {
    // target is mapping query, which needs to exist to create psl.
    if (not fTransMaps[0]->haveQuerySeq(features[0]->getSeqid())) {
        return NULL;
    }
    int tSize = fTransMaps[0]->getQuerySeqSize(features[0]->getSeqid()); // target is mapping query
    struct psl* srcPsl = FeaturesToPsl::toPsl(qName, tSize, features);
    PslVector mappedPsls = recursiveMapPsl(srcPsl, 0);
    return new PslMapping(srcPsl, mappedPsls);
}


/* map a single feature */
PslMapping* FeatureTransMap::mapFeature(const string& qName,
                                        const Feature* feature) const {
    FeatureVector features;
    features.push_back(const_cast<Feature*>(feature));
    return mapFeatures(qName, features);
}
