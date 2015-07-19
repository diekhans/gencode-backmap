
#include "geneMapper.hh"
#include "jkcommon.hh"
extern "C" {
#define hash jkhash
#define new jknew
#include "psl.h"
#include "dnautil.h"
#undef hash
#undef new
}
#include "transMapper.hh"
#include "gxf.hh"

/** compute total size of features */
int GeneMapper::sumFeatureSizes(const GxfFeatureVector& features) {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->fEnd - features[i]->fStart)+1;
    }
    return sz;
}

/* create a psl from a list of features. assumes features are sorter in
 * ascending order */
struct psl* GeneMapper::featuresToPsl(const string& qName,
                                      const GxfFeatureVector& exons) {
    // this does a [1..n] to [0..n) conversion
    int qSize = sumFeatureSizes(exons);
    int tSize = fTransMapper->getQuerySeqSize(exons[0]->fSeqid); // target is mapping query
    struct psl* psl = pslNew(toCharStr(qName), qSize, 0, qSize,
                             toCharStr(exons[0]->fSeqid), tSize, exons[0]->fStart-1, exons[exons.size()-1]->fEnd,
                             toCharStr(exons[0]->fStrand), exons.size(), 0);
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < exons.size(); iBlk++) {
        const GxfFeature* feature = exons[iBlk];
        psl->blockSizes[iBlk] = (feature->fEnd - feature->fStart)+1;
        psl->qStarts[iBlk] = qStart;
        psl->tStarts[iBlk] = feature->fStart-1;
        psl->match += psl->blockSizes[iBlk];
        qStart += psl->blockSizes[iBlk];
        if ((iBlk > 0) && (psl->tStarts[iBlk] > (psl->tStarts[iBlk-1]+psl->blockSizes[iBlk]))) {
            throw invalid_argument("out of order blocks: " + feature->toString());
        }
    }
    return psl;
}

/* get exon features */
GxfFeatureVector GeneMapper::getExons(const GxfFeatureNode* transcript) {
    GxfFeatureVector exons;
    for (size_t i = 0; i < transcript->fChildren.size(); i++) {
        if (transcript->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            exons.push_back(transcript->fChildren[i]->fFeature);
        }
    }
    return exons;
}

/* create a psl for exons of a transcript */
struct psl* GeneMapper::transcriptExonsToPsl(const GxfFeatureNode* transcript) {
    const string& qName = transcript->fFeature->getAttr("transcript_id")->fVal;
    const GxfFeatureVector exons = getExons(transcript);
    return featuresToPsl(qName, exons);
        
}
