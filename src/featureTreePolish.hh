/*
 * Post mapping corrections to feature tree.
 */
#ifndef featureTreePolish_hh
#define featureTreePolish_hh
#include <assert.h>
#include <map>
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
    // use to map old exon ids to new exons in a transcript
    typedef map<int, vector<FeatureNode*> > ExonIdExonMap;
    typedef ExonIdExonMap::iterator ExonIdExonMapIter;
    typedef ExonIdExonMap::const_iterator ExonIdExonMapConstIter;

    static void renumberExon(FeatureNode* exonNode,
                             int exonNum,
                             ExonIdExonMap& exonIdExonMap);
    static void renumberExons(FeatureNode* transcriptTree,
                              ExonIdExonMap& exonIdExonMap);
    static FeatureNode* findNewExon(FeatureNode* featureNode,
                                    int oldExonNum,
                                    ExonIdExonMap& exonIdExonMap);
    static void renumberOtherFeature(FeatureNode* featureNode,
                                     ExonIdExonMap& exonIdExonMap);
    static void renumberOtherFeatures(FeatureNode* featureNode,
                                      ExonIdExonMap& exonIdExonMap);
    static void renumberTranscript(FeatureNode* transcriptTree);
    static void renumberGeneExons(FeatureNode* geneTreeRoot);

    public:
    /* last minute fix-ups */
    static void polishGene(FeatureNode* geneTreeRoot);
};

#endif
