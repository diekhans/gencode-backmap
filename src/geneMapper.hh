/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "gxfFeatureTree.hh"
class TransMap;
class PslMapping;
struct psl;
class PslCursor;

/* class that maps a gene to the new assemble */
class GeneMapper {
    private:
    const TransMap* fGenomeTransMap;
    
    void processGene(GxfParser *gxfParser,
                     const GxfFeature* geneFeature,
                     ostream& outFh) const;
    void processTranscript(const GxfFeatureNode* transcriptTree) const;

    public:
    /* Constructor */
    GeneMapper(const TransMap* genomeTransMap):
        fGenomeTransMap(genomeTransMap) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfParser *gxfParser,
                ostream& outFh) const;
};

#endif
