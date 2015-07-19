/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
class TransMapper;
struct psl;

/* class that maps a gene to the new assemble */
class GeneMapper {
    private:
    const TransMapper* fTransMapper;  // object to performance mappings
    
    int sumFeatureSizes(const GxfFeatureVector& features);
    struct psl* featuresToPsl(const string& qName,
                              const GxfFeatureVector& exons);
    GxfFeatureVector getExons(const GxfFeatureNode* transcript);
    struct psl* transcriptExonsToPsl(const GxfFeatureNode* transcript);
};

#endif
