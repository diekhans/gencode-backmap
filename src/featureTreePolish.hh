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

    void renumberExon(FeatureNode* exonNode,
                      int exonNum,
                      ExonIdExonMap& exonIdExonMap) const;
    void renumberExons(FeatureNode* transcriptTree,
                       ExonIdExonMap& exonIdExonMap) const;
    FeatureNode* findNewExon(FeatureNode* featureNode,
                             int oldExonNum,
                             ExonIdExonMap& exonIdExonMap) const;
    void renumberOtherFeature(FeatureNode* featureNode,
                              ExonIdExonMap& exonIdExonMap) const;
    void renumberOtherFeatures(FeatureNode* featureNode,
                               ExonIdExonMap& exonIdExonMap) const;
    void renumberTranscriptExons(FeatureNode* transcriptTree) const;
    void renumberGeneExons(FeatureNode* geneTree) const;
    bool isRemapped(FeatureNode* featureNode) const;
    void setMappingVersionInId(FeatureNode* featureNode,
                               const AttrVal* attr,
                               int version) const;
    void setMappingVersion(FeatureNode* featureNode,
                           const string& idAttrName,
                           const string& havanaIdAttrName,
                           int version) const;
    void recursiveSetMappingVersion(FeatureNode* featureNode,
                                    const string& idAttrName,
                                    const string& havanaIdAttrName,
                                    int version) const;
    void setTranscriptMappingVersion(FeatureNode* transcriptTree) const;
    void setGeneMappingVersion(FeatureNode* geneTree) const;
    
    public:
    /* last minute fix-ups */
    void polishGene(FeatureNode* geneTree) const;
};

#endif
