/*
 * Container for a source and resulting mapped PSLs.
 */
#ifndef pslMapping_hh
#define pslMapping_hh
#include "pslOps.hh"

/* set of mapped alignments */
class PslMapping {
    public:
    struct psl* fSrcPsl;
    PslVector fMappedPsls;  // top one is best scoring
    int fScore;  // mapping score of top psl; 0 is perfect

    private:
    static int numPslAligned(const struct psl* psl);
    void sortMappedPsls();

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls);

    /* free up psls */
    ~PslMapping();

    /* Compute a mapping score between the src and mapped psl.  A perfect
     * mapping is a zero score.  Extra inserts count against the score. */
    static int calcPslMappingScore(const struct psl* srcPsl,
                                   const struct psl* mappedPsl);
};


#endif
