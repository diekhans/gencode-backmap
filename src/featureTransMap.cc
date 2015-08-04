/*
 * Mapping of GxF features using transmap
 */
#include "featureTransMap.hh"
#include "gxf.hh"

/** compute total size of features */
int FeatureTransMap::sumFeatureSizes(const GxfFeatureVector& features) const {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->fEnd - features[i]->fStart)+1;
    }
    return sz;
}

/* build blocks for positive strand alignment. */
void FeatureTransMap::makePosPslBlocks(struct psl* psl,
                                  const GxfFeatureVector& exons) const {
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < exons.size(); iBlk++) {
        const GxfFeature* feature = exons[iBlk];
        psl->blockSizes[iBlk] = (feature->fEnd - feature->fStart)+1;
        psl->qStarts[iBlk] = qStart;
        psl->tStarts[iBlk] = feature->fStart-1;
        psl->match += psl->blockSizes[iBlk];
        qStart += psl->blockSizes[iBlk];
        if ((iBlk > 0) && (psl->tStarts[iBlk] < (psl->tStarts[iBlk-1]+psl->blockSizes[iBlk-1]))) {
            throw invalid_argument("out of order blocks: " + feature->toString());
        }
        psl->blockCount++;
    }
}

/* build blocks for negative strand alignment. */
void FeatureTransMap::makeNegPslBlocks(struct psl* psl,
                                  const GxfFeatureVector& exons) const {
    int qStart = 0;
    size_t iExon = exons.size() - 1;
    for (size_t iBlk = 0; iBlk < exons.size(); iBlk++, iExon--) {
        const GxfFeature* feature = exons[iExon];
        psl->blockSizes[iBlk] = (feature->fEnd - feature->fStart)+1;
        psl->qStarts[iBlk] = psl->qSize - qStart;
        psl->tStarts[iBlk] = feature->fStart-1;
        psl->match += psl->blockSizes[iBlk];
        qStart += psl->blockSizes[iBlk];
        if ((iBlk > 0) && (psl->tStarts[iBlk] < (psl->tStarts[iBlk-1]+psl->blockSizes[iBlk-1]))) {
            throw invalid_argument("out of order blocks: " + feature->toString());
        }
        psl->blockCount++;
    }
}

/* create a psl from a list of features. assumes features are sorter in
 * ascending order */
struct psl* FeatureTransMap::featuresToPsl(const string& qName,
                                           const GxfFeatureVector& features) const {
    // this does a [1..n] to [0..n) conversion
    int qSize = sumFeatureSizes(features);
    int tSize = fTransMap->getQuerySeqSize(features[0]->fSeqid); // target is mapping query
    struct psl* psl = pslNew(toCharStr(qName), qSize, 0, qSize,
                             toCharStr(features[0]->fSeqid), tSize, features[0]->fStart-1, features[features.size()-1]->fEnd,
                             toCharStr(features[0]->fStrand), features.size(), 0);
    if (psl->strand[0] == '+') {
        makePosPslBlocks(psl, features);
    } else {
        makeNegPslBlocks(psl, features);
    }
    return psl;
}

/* map features to list of PSLs */
PslMapping* FeatureTransMap::mapFeatures(const string& qName,
                                         const GxfFeatureVector& features) const {
    struct psl* inPsl = featuresToPsl(qName, features);
    return fTransMap->mapPsl(inPsl);
}
