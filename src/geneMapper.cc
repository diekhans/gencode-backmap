
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
#include "gxf.hh"
#include <algorithm>
#include <fstream>
#include <iostream>

/* constructor, sort mapped PSLs */
PslMapping::PslMapping(struct psl* srcPsl,
                       PslVector& mappedPsls):
    fSrcPsl(srcPsl),
    fMappedPsls(mappedPsls),
    fScore(0) {
    sortMappedPsls();
    fScore = calcPslMappingScore(srcPsl, mappedPsls[0]); // best score
}

/* free up psls */
PslMapping::~PslMapping() {
    pslFree(&fSrcPsl);
    for (size_t i = 0; i < fMappedPsls.size(); i++) {
        pslFree(&(fMappedPsls[i]));
    }
}

/* calculate the number of aligned bases */
int PslMapping::numPslAligned(struct psl* psl) {
    return (psl->match + psl->repMatch + psl->misMatch);
}

/* Compute a mapping score between a src and mapped psl.  A perfect mapping is
 * a zero score.  Extra inserts count against the score.
 */
int PslMapping::calcPslMappingScore(struct psl* srcPsl,
                                    struct psl* mappedPsl) {
    return (numPslAligned(srcPsl) - numPslAligned(mappedPsl))
        + abs(int(srcPsl->qNumInsert)-int(mappedPsl->qNumInsert))
        + abs(int(srcPsl->tNumInsert)-int(mappedPsl->tNumInsert));
}

/* comparison functor based on lowest store. */
class ScoreCmp {
    public:
    struct psl* fSrcPsl;
    
    ScoreCmp(struct psl* srcPsl):
        fSrcPsl(srcPsl) {
    }

    int operator()(struct psl* mappedPsl1,
                   struct psl* mappedPsl2) const {
        return -(PslMapping::calcPslMappingScore(fSrcPsl, mappedPsl1)
                 - PslMapping::calcPslMappingScore(fSrcPsl, mappedPsl2));
    }
};
/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls() {
    ScoreCmp scoreCmp(fSrcPsl);
    sort(fMappedPsls.begin(), fMappedPsls.end(), scoreCmp);
}


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
    int tSize = fTransMap->getQuerySeqSize(exons[0]->fSeqid); // target is mapping query
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
GxfFeatureVector GeneMapper::getExons(const GxfFeatureNode* transcriptTree) {
    GxfFeatureVector exons;
    for (size_t i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            exons.push_back(transcriptTree->fChildren[i]->fFeature);
        }
    }
    return exons;
}

/* create a psl for exons of a transcript */
struct psl* GeneMapper::transcriptExonsToPsl(const GxfFeatureNode* transcriptTree) {
    const string& qName = transcriptTree->fFeature->getAttr("transcript_id")->fVal;
    const GxfFeatureVector exons = getExons(transcriptTree);
    return featuresToPsl(qName, exons);
}

/* 
 * Map a transcript's exons, returning object that has mappings and scoring
 */
PslMapping* GeneMapper::mapTranscriptExons(const GxfFeatureNode* transcriptTree) {
    struct psl* transcriptPsl = transcriptExonsToPsl(transcriptTree);
    PslVector mappedTranscriptPsls = fTransMap->mapPsl(transcriptPsl);
    return new PslMapping(transcriptPsl, mappedTranscriptPsls);
}


/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfParser *gxfParser,
                        ostream& outFh) {
    const GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (instanceOf(gxfRecord, GxfFeature)) {
            const GxfFeature* geneFeature = dynamic_cast<const GxfFeature*>(gxfRecord);
            GxfFeatureTree* geneTree = new GxfFeatureTree(gxfParser, geneFeature);
            geneTree->write(outFh);
            delete geneTree;
        } else {
            outFh << gxfRecord->toString() << endl;
            delete gxfRecord;
        }
    }
}
