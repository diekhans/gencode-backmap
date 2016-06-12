/*
 * Container for a source and resulting mapped PSLs.
 */
#include "pslMapping.hh"
#include <algorithm>
#include "gxfRecord.hh"
#include <iostream>

// FIXME: passing down features to this level in simple container is annoying.
// It would be better to have a sort function passed it, but we had all the
// core dump issues with sort, put algorithm here for now.

/* constructor, sort mapped PSLs */
PslMapping::PslMapping(struct psl* srcPsl,
                       PslVector& mappedPsls,
                       const GxfFeature* primaryTarget,
                       const GxfFeature* secondaryTarget):
    fSrcPsl(srcPsl),
    fMappedPsl(NULL),
    fMappedPsls(mappedPsls) {
    assert(pslQStrand(srcPsl) == '+');
    if (fMappedPsls.size() > 0) {
        sortMappedPsls(primaryTarget, secondaryTarget);
    }
}

/* free up psls */
PslMapping::~PslMapping() {
    pslFree(&fSrcPsl);
    for (size_t i = 0; i < fMappedPsls.size(); i++) {
        pslFree(&(fMappedPsls[i]));
    }
}

/* calculate the number of aligned bases */
int PslMapping::numAlignedBases(const struct psl* psl) {
    return (psl->match + psl->repMatch + psl->misMatch);
}

/* Compute a mapping score between a src and mapped psl.  A perfect mapping is
 * a zero score.  Extra inserts count against the score.
 */
int PslMapping::calcPslMappingScore(const struct psl* srcPsl,
                                    const struct psl* mappedPsl) {
    // FIXME: not sure if this makes sens
    return (numAlignedBases(srcPsl) - numAlignedBases(mappedPsl))
        + abs(int(srcPsl->qNumInsert)-int(mappedPsl->qNumInsert))
        + abs(int(srcPsl->tNumInsert)-int(mappedPsl->tNumInsert));
}

/* dump for debugging purposes, adding optional description */
void PslMapping::dump(ostream& fh,
                      const string& description,
                      const string& indent) const {
    if (description != "") {
        fh  << description << endl;
    }
    fh << indent << "srcPsl:       " << pslToString(fSrcPsl) << endl;
    if (fMappedPsls.size() == 0) {
        fh << indent << "mappedPsl[0]: none" << endl;
    } else {
        for (int i = 0; i < fMappedPsls.size(); i++) {
            fh << indent << "mappedPsl[" << i << "]: " << pslToString(fMappedPsls[i]) << endl;
        }
    }
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


#if 0
// FIXME: with this seemingly correct sort, the C++ library goes of the end of the vector and SEGVs
// this happens on g++ 4.8.2 on Linux and 4.9 on OS/X (Mac ports).  It doesn't happen with
// clang on OS/X.
/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls() {
    //sort(fMappedPsls.begin(), fMappedPsls.end(), ScoreCmp(fSrcPsl));
    std::sort(fMappedPsls.begin(), fMappedPsls.end(),
              [this](const struct psl* a, const struct psl* b) -> bool {
                  return (PslMapping::calcPslMappingScore(this->fSrcPsl, a)
                          - PslMapping::calcPslMappingScore(this->fSrcPsl, b));
              });
    if (fMappedPsls.size() > 0) {
        fMappedPsl = fMappedPsls[0];
    }
}
#else
/* globals for use in comparison because qsort doesn't have a client data */
static struct psl* gSrcPsl = NULL;
static const GxfFeature* gPrimaryTarget = NULL;
static const GxfFeature* gSecondaryTarget = NULL;

/* compute fraction of overlap similarity for a psl and a target feature. */
static float targetSimilarity(const struct psl *mappedPsl,
                              const GxfFeature* targetFeature) {
    if (targetFeature->getSeqid() != mappedPsl->tName) {
        return 0.0;  // different chrom
    }
    int maxStart = max(targetFeature->getStart()-1, mappedPsl->tStart);
    int minEnd = max(targetFeature->getEnd(), mappedPsl->tEnd);
    if (minEnd <= maxStart) {
        return 0.0;  // no overlap
    }
    return float(2*(minEnd - maxStart)) / float((mappedPsl->tEnd-mappedPsl->tStart) + targetFeature->size());
}

/* compare function based on target similarity. */
static int targetSimilarityCmp(const struct psl *mappedPsl1,
                               const struct psl *mappedPsl2,
                               const GxfFeature* targetFeature) {
    // don't think we need an approximate compare, because it will be 0.0 if no overlap,
    // and don't know why very close overlap would happen
    float sim1 = targetSimilarity(mappedPsl1, targetFeature);
    float sim2 = targetSimilarity(mappedPsl2, targetFeature);
    if (sim1 < sim2) {
        return -1;
    } else if (sim1 > sim2) {
        return 1;
    } else {
        return 0;
    }
}

/* compare by span similarity */
static int spanSimilarityCmp(const struct psl *mappedPsl1,
                             const struct psl *mappedPsl2) {
    int srcSpan = (gSrcPsl->tEnd - gSrcPsl->tStart);
    int span1Diff = abs(srcSpan - (mappedPsl1->tEnd - mappedPsl1->tStart));
    int span2Diff = abs(srcSpan - (mappedPsl2->tEnd - mappedPsl2->tStart));
    // rank smallest span change best (reverse sort)
    if (span1Diff < span2Diff) {
        return -1;
    } else if (span1Diff > span2Diff) {
        return 1;
    } else {
        return 0;
    }
}

/* compare two psl to see which is better mapped. */
static int mapScoreCmp(const void *va, const void *vb) {
    const struct psl *mappedPsl1 = *static_cast<const struct psl *const*>(va);
    const struct psl *mappedPsl2 = *static_cast<const struct psl *const*>(vb);
    if (gPrimaryTarget != NULL) {
        int diff = targetSimilarityCmp(mappedPsl1, mappedPsl2, gPrimaryTarget);
        if (diff != 0) {
            return -diff; // inverse sort
        }
    }
    if (gSecondaryTarget != NULL) {
        int diff = targetSimilarityCmp(mappedPsl1, mappedPsl2, gSecondaryTarget);
        if (diff != 0) {
            return -diff; // inverse sort
        }
    }
    int diff = spanSimilarityCmp(mappedPsl1, mappedPsl2);
    if (diff != 0) {
        return diff;
    }
    // FIXME: this maybe silly
    return (PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl1)
            - PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl2));
}

/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls(const GxfFeature* primaryTarget,
                                const GxfFeature* secondaryTarget) {
    gSrcPsl = fSrcPsl;
    gPrimaryTarget = primaryTarget;
    gSecondaryTarget = secondaryTarget;
    qsort(static_cast<struct psl**>(&(fMappedPsls[0])), fMappedPsls.size(), sizeof(struct psl*), mapScoreCmp);
    gSrcPsl = NULL;
    gPrimaryTarget = gSecondaryTarget = NULL;
    if (fMappedPsls.size() > 0) {
        fMappedPsl = fMappedPsls[0];
    }
}
#endif



