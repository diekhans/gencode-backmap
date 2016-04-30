/*
 * Post mapping corrections to feature tree.
 */
#ifndef featureTreePolish_hh
#define featureTreePolish_hh
#include <assert.h>
#include "featureTree.hh"

/*
 * Notes:
 *  - this class will also handle gap filling in the future.
 */

/*
 * Post mapping corrections to feature tree.
 */
class FeatureTreePolish {
    private:

    static void renumberExon(GxfFeature* exon,
                             int exonNum);
    static void renumberExons(FeatureNode* transcriptTree);
    static void renumberGeneExons(FeatureNode* geneTreeRoot);

    public:
    /* last minute fix-ups */
    static void polishGene(FeatureNode* geneTreeRoot);
};

#endif
