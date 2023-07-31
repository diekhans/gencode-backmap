/*
 * Container for a source and resulting mapped PSLs.
 */
#include "pslMapping.hh"
#include <algorithm>
#include "featureTree.hh"
#include <iostream>

// FIXME: passing down features to this level in simple container is annoying.
// It would be better to have a sort function passed it, but we had all the
// core dump issues with sort, put algorithm here for now.
// NOTE: core dump was due to not following ordering requirements of C++ sort.

static const bool debug = false;


#undef USE_CXX_SORT

/* constructor, sort mapped PSLs */
PslMapping::PslMapping(struct psl* srcPsl,
                       PslVector& mappedPsls):
    fSrcPsl(srcPsl),
    fMappedPsls(mappedPsls) {
    assert(pslQStrand(srcPsl) == '+');
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
    // this mostly sort out ties, not that meaningful
    return abs(numAlignedBases(srcPsl) - numAlignedBases(mappedPsl))
        + abs(int(srcPsl->qNumInsert) - int(mappedPsl->qNumInsert))
        + abs(int(srcPsl->tNumInsert) - int(mappedPsl->tNumInsert));
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

static void dumpMapped(struct psl* srcPsl,
                       PslVector& mappedPsls) {
#if debug
    if (mappedPsls.size() > 1) {
        fprintf(stdout, "@ ");
        pslTabOut(srcPsl, stdout);
        for (auto it: mappedPsls) {
            pslTabOut(it, stdout);
        }
        fflush(stdout);
    }
#endif
}

/* compute fraction of overlap similarity for a psl and a target feature. */
static float targetSimilarity(const struct psl *mappedPsl,
                              const FeatureNode* targetFeature) {
    if (targetFeature->getSeqid() != mappedPsl->tName) {
        return 0.0;  // different chrom
    }
    int maxStart = max(targetFeature->getStart()-1, mappedPsl->tStart);
    int minEnd = max(targetFeature->getEnd(), mappedPsl->tEnd);
    if (minEnd <= maxStart) {
        return 0.0;  // no overlap
    }
    return float(2*(minEnd - maxStart)) / float((mappedPsl->tEnd-mappedPsl->tStart) + targetFeature->length());
}

/* compare function based on target similarity, higher is better (-1)  */
static int targetSimilarityCmp(const struct psl *mappedPsl1,
                               const struct psl *mappedPsl2,
                               const FeatureNode* targetFeature) {
    // don't think we need an approximate compare, because it will be 0.0 if no overlap,
    // and don't know why very close overlap would happen
    float sim1 = targetSimilarity(mappedPsl1, targetFeature);
    float sim2 = targetSimilarity(mappedPsl2, targetFeature);
    if (sim1 < sim2) {
        return 1;  // better
    } else if (sim1 > sim2) {
        return -1;
    } else {
        return 0;
    }
}

/* compare by span similarity */
static int spanSimilarityCmp(const struct psl *mappedPsl1,
                             const struct psl *mappedPsl2,
                             const struct psl *srcPsl) {
    int srcSpan = (srcPsl->tEnd - srcPsl->tStart);
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

#ifdef USE_CXX_SORT
/* std::sort compare implementation */
static bool mappedPslCompare(struct psl* srcPsl,
                             const FeatureNode* primaryTarget,
                             const FeatureNode* secondaryTarget,
                             const struct psl* mappedPsl1,
                             const struct psl* mappedPsl2) {
    // must provide strict less than
    if (mappedPsl1 == mappedPsl2) {
        return false;  // same entry, hope this avoids corruption on bad order
    }
    // reverse sort, less than is better
    if ((primaryTarget != NULL)
        and (targetSimilarityCmp(mappedPsl1, mappedPsl2, primaryTarget) < 0)) {
        return true;
    }
    if ((secondaryTarget != NULL)
        and (targetSimilarityCmp(mappedPsl1, mappedPsl2, secondaryTarget) < 0)) {
        return true;
    }
    if (spanSimilarityCmp(mappedPsl1, mappedPsl2, srcPsl) < 0) {
        return true;
    }
    if (PslMapping::calcPslMappingScore(srcPsl, mappedPsl1) <
        PslMapping::calcPslMappingScore(srcPsl, mappedPsl2)) {
        return true;
    }
    return false;
}

/* check it strict weak for a pair */
static void checkStrictWeak(struct psl* srcPsl,
                            const FeatureNode* primaryTarget,
                            const FeatureNode* secondaryTarget,
                            const struct psl* mappedPsl1,
                            const struct psl* mappedPsl2) {
    if (mappedPslCompare(srcPsl, primaryTarget, secondaryTarget, mappedPsl1, mappedPsl1)) {
        throw logic_error("PSL1 less than self: " + pslToString(mappedPsl1));
    }
    if (mappedPslCompare(srcPsl, primaryTarget, secondaryTarget, mappedPsl2, mappedPsl2)) {
        throw logic_error("PSL2 less than self: " + pslToString(mappedPsl2));
    }
    bool xyLessThan = mappedPslCompare(srcPsl, primaryTarget, secondaryTarget, mappedPsl1, mappedPsl2);
    bool yxLessThan = mappedPslCompare(srcPsl, primaryTarget, secondaryTarget, mappedPsl2, mappedPsl1);
    if (xyLessThan and yxLessThan) {
        throw logic_error("PSL1/PSL2 asymmetry error: " + pslToString(mappedPsl1)
                          + " and " + pslToString(mappedPsl2));
    }
}

/*
 * Validate that strict weak ordering holds:
 *   - For all x in S, it is not the case that x < x (irreflexivity).
 *   - For all x, y in S, if x < y then it is not the case that y < x (asymmetry).
 *   - For all x, y, z in S, if x < y and y < z then x < z (transitivity).
 *   - For all x, y, z in S, if x is incomparable with y (neither x < y
 *     nor y < x hold), and y is incomparable with z, then x is
 *     incomparable with z (transitivity of incomparability).
 * only the first two are check
 */
static void checkStrictWeak(struct psl* srcPsl,
                            const FeatureNode* primaryTarget,
                            const FeatureNode* secondaryTarget,
                            const PslVector& mappedPsls) {
    for (auto it1 = mappedPsls.begin(); it1 != mappedPsls.end(); it1++) {
        for (auto it2 = it1 + 1; it2 != mappedPsls.end(); it2++) {
            checkStrictWeak(srcPsl,primaryTarget, secondaryTarget, *it1, *it2);
        }
    }
}

/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls(const FeatureNode* primaryTarget,
                                const FeatureNode* secondaryTarget) {
    dumpMapped(fSrcPsl, fMappedPsls);
    checkStrictWeak(fSrcPsl, primaryTarget, secondaryTarget, fMappedPsls);
    std::sort(fMappedPsls.begin(), fMappedPsls.end(),
              [this, primaryTarget, secondaryTarget](const struct psl* mappedPsl1, const struct psl* mappedPsl2) -> bool {
                  return mappedPslCompare(this->fSrcPsl, primaryTarget, secondaryTarget,
                                          mappedPsl1, mappedPsl2);
              });
    if (fMappedPsls.size() > 0) {
        fMappedPsl = fMappedPsls[0];
    }
}
#else
/* globals for use in comparison because qsort doesn't have a client data */
static struct psl* gSrcPsl = NULL;
static const FeatureNode* gPrimaryTarget = NULL;
static const FeatureNode* gSecondaryTarget = NULL;

/* compare two psl to see which is better mapped. */
static int mapScoreCmp(const void *va, const void *vb) {
    const struct psl *mappedPsl1 = *static_cast<const struct psl *const*>(va);
    const struct psl *mappedPsl2 = *static_cast<const struct psl *const*>(vb);
    if (gPrimaryTarget != NULL) {
        int diff = targetSimilarityCmp(mappedPsl1, mappedPsl2, gPrimaryTarget);
        if (diff != 0) {
            return diff;
        }
    }
    if (gSecondaryTarget != NULL) {
        int diff = targetSimilarityCmp(mappedPsl1, mappedPsl2, gSecondaryTarget);
        if (diff != 0) {
            return diff;
        }
    }
    int diff = spanSimilarityCmp(mappedPsl1, mappedPsl2, gSrcPsl);
    if (diff != 0) {
        return diff;
    }
    // FIXME: this maybe silly
    return (PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl1)
            - PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl2));
}

/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls(const FeatureNode* primaryTarget,
                                const FeatureNode* secondaryTarget) {
    dumpMapped(fSrcPsl, fMappedPsls);
    gSrcPsl = fSrcPsl;
    gPrimaryTarget = primaryTarget;
    gSecondaryTarget = secondaryTarget;
    qsort(static_cast<struct psl**>(&(fMappedPsls[0])), fMappedPsls.size(), sizeof(struct psl*), mapScoreCmp);
    gSrcPsl = NULL;
    gPrimaryTarget = gSecondaryTarget = NULL;
}
#endif

/* filter for mappings to same mapped target sequence as source */
void PslMapping::filterSameTarget() {
    for (auto iter = fMappedPsls.begin(); iter != fMappedPsls.end();) {
        struct psl *mpsl = *iter;
        if (!sameString(mpsl->tName, fSrcPsl->tName)) {
            pslFree(&mpsl);
            iter = fMappedPsls.erase(iter);
        } else {
            iter++;
        }
    }
}


