/*
 * Container for a source and resulting mapped PSLs.
 */
#include "pslMapping.hh"
#include <algorithm>

/* constructor, sort mapped PSLs */
PslMapping::PslMapping(struct psl* srcPsl,
                       PslVector& mappedPsls):
    fSrcPsl(srcPsl),
    fMappedPsls(mappedPsls),
    fScore(-1) {
    if (fMappedPsls.size() > 0) {
        sortMappedPsls();
        fScore = calcPslMappingScore(srcPsl, fMappedPsls[0]); // best score due to sort
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
int PslMapping::numPslAligned(const struct psl* psl) {
    return (psl->match + psl->repMatch + psl->misMatch);
}

/* Compute a mapping score between a src and mapped psl.  A perfect mapping is
 * a zero score.  Extra inserts count against the score.
 */
int PslMapping::calcPslMappingScore(const struct psl* srcPsl,
                                    const struct psl* mappedPsl) {
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


#if 0
// FIXME: with this seemingly correct sort, the C++ library goes of the end of the vector and SEGVs
// this happens on g++ 4.8.2 on Linux and 4.9 on OS/X (Mac ports).  It doesn't happen with
// clang on OS/X.
/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls() {
    //sort(fMappedPsls.begin(), fMappedPsls.end(), ScoreCmp(fSrcPsl));
    std::sort(fMappedPsls.begin(), fMappedPsls.end(),
              [this](const struct psl* a, const struct psl* b) -> bool {
                  return -(PslMapping::calcPslMappingScore(this->fSrcPsl, a)
                           - PslMapping::calcPslMappingScore(this->fSrcPsl, b));
              });
}
#else
/* compare two psl to see which is better mapped. */
static struct psl* gSrcPsl = NULL;
static int mapScoreCmp(const void *va, const void *vb) {
    const struct psl *mappedPsl1 = *static_cast<const struct psl *const*>(va);
    const struct psl *mappedPsl2 = *static_cast<const struct psl *const*>(vb);
    return -(PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl1)
             - PslMapping::calcPslMappingScore(gSrcPsl, mappedPsl2));
}

/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls() {
    gSrcPsl = fSrcPsl;
    qsort(static_cast<struct psl**>(&(fMappedPsls[0])), fMappedPsls.size(), sizeof(struct psl*), mapScoreCmp);
    gSrcPsl = NULL;
}
#endif



