/*
 * Container for a source and resulting mapped PSLs.
 */
#ifndef pslMapping_hh
#define pslMapping_hh
#include "pslOps.hh"
#include <iostream>
class GxfFeature;

/* set of mapped alignments */
class PslMapping {
    private:
    struct psl* fSrcPsl;
    struct psl* fMappedPsl; // best mapped PSL, or NULL if none mapped
    PslVector fMappedPsls;  // all mapped PSLs, [0] is fMappedPsl

    static int numAlignedBases(const struct psl* psl);

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls,
               const GxfFeature* primaryTarget=NULL,
               const GxfFeature* secondaryTarget=NULL);

    /* free up psls */
    ~PslMapping();

    /* Accessors */
    struct psl* getSrcPsl() const {
        return fSrcPsl;
    }
    struct psl* getMappedPsl() const {
        return fMappedPsl;
    }
    const PslVector getMappedPsls() const {
        return fMappedPsls;
    }

    /* sort or restore the mapped PSLs  Kind of a hack to allow restore so
     * we don't have to pass primaryTarget/secondaryTargets everywhere */
    void sortMappedPsls(const GxfFeature* primaryTarget=NULL,
                        const GxfFeature* secondaryTarget=NULL);

    /* are there any mappings? */
    bool haveMappings() const {
        return fMappedPsls.size() > 0;
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
