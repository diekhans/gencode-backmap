/*
 * Container for a source and resulting mapped PSLs.
 */
#ifndef pslMapping_hh
#define pslMapping_hh
#include "pslOps.hh"
#include <iostream>
class FeatureNode;

/* set of mapped alignments */
class PslMapping {
    private:
    struct psl* fSrcPsl;
    PslVector fMappedPsls;  // all mapped PSLs, [0] is best mapped after sort

    static int numAlignedBases(const struct psl* psl);

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls);

    /* free up psls */
    ~PslMapping();

    /* Accessors */
    struct psl* getSrcPsl() const {
        return fSrcPsl;
    }
    struct psl* getBestMappedPsl() const {
        return fMappedPsls.empty() ? NULL: fMappedPsls[0];
    }
    const PslVector getMappedPsls() const {
        return fMappedPsls;
    }

    /* Sort the mapped PSLs by score, optionally using primaryTarget/secondaryTargets */
    void sortMappedPsls(const FeatureNode* primaryTarget=NULL,
                        const FeatureNode* secondaryTarget=NULL);

    /* filter for mappings to same mapped target sequence as source */
    void filterSameTarget();

    /* are there any mappings? */
    bool haveMappings() const {
        return not fMappedPsls.empty();
    }

    /** write the mapped PSLs */
    void writeMapped(ostream& fh) const {
        for (size_t i = 0; i < fMappedPsls.size(); i++) {
            fh << pslToString(fMappedPsls[i]) << endl;
        }
    }
    
    /* dump for debugging purposes, adding optional description */
    void dump(ostream& fh,
              const string& description="",
              const string& indent="") const;
    
    /* Compute a mapping score between the src and mapped psl.  A perfect
     * mapping is a zero score.  Extra inserts count against the score. */
    static int calcPslMappingScore(const struct psl* srcPsl,
                                   const struct psl* mappedPsl);
};


#endif
