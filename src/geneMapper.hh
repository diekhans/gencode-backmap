/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "transMapper.hh"
#include "gxfFeatureTree.hh"
struct psl;

/* set of mapped alignments */
class PslMapping {
    public:
    struct psl* fSrcPsl;
    PslVector fMappedPsls;  // top one is best scoring
    int fScore;  // mapping score of top psl; 0 is perfect

    private:
    static int numPslAligned(struct psl* psl);
    void sortMappedPsls();

    public:
    /* constructor, sort mapped PSLs */
    PslMapping(struct psl* srcPsl,
               PslVector& mappedPsls);

    /* free up psls */
    ~PslMapping();

    /* Compute a mapping score between the src and mapped psl.  A perfect
     * mapping is a zero score.  Extra inserts count against the score. */
    static int calcPslMappingScore(struct psl* srcPsl,
                                   struct psl* mappedPsl);
};

/* class that maps a gene to the new assemble */
class GeneMapper {
    private:
    const TransMapper* fTransMapper;  // object to performance mappings
    
    int sumFeatureSizes(const GxfFeatureVector& features);
    struct psl* featuresToPsl(const string& qName,
                              const GxfFeatureVector& exons);
    GxfFeatureVector getExons(const GxfFeatureNode* transcript);
    struct psl* transcriptExonsToPsl(const GxfFeatureNode* transcriptTree);
    PslMapping* mapTranscriptExons(const GxfFeatureNode* transcriptTree);

    public:
    /* Constructor */
    GeneMapper(const TransMapper* transMapper):
        fTransMapper(transMapper) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfParser *gxfParser,
                ostream& outFh);
};

#endif
