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
    // use to map old exon numbers to new exons in a transcript
    typedef map<int, vector<FeatureNode*> > ExonNumExonMap;
    typedef ExonNumExonMap::iterator ExonNumExonMapIter;
    typedef ExonNumExonMap::const_iterator ExonNumExonMapConstIter;

    // use to add mapping version number to exon ids
    typedef map<string, vector<FeatureNode*> > ExonIdExonMap;
    typedef ExonIdExonMap::iterator ExonIdExonMapIter;
    typedef ExonIdExonMap::const_iterator ExonIdExonMapConstIter;
    
    void renumberExon(FeatureNode* exonNode,
                      int exonNum,
                      ExonNumExonMap& exonNumExonMap) const;
    void renumberExons(FeatureNode* transcriptTree,
                       ExonNumExonMap& exonNumExonMap) const;
    FeatureNode* findNewExon(FeatureNode* featureNode,
                             int oldExonNum,
                             ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeature(FeatureNode* featureNode,
                              ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeatures(FeatureNode* featureNode,
                               ExonNumExonMap& exonNumExonMap) const;
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
    void recordTranscriptMappedExons(FeatureNode* transcriptTree,
                                     ExonIdExonMap& exonIdExonMap) const;
    void setTranscriptMappingVersion(FeatureNode* transcriptTree,
                                     ExonIdExonMap& exonIdExonMap) const;
    void setExonMappingVersion(FeatureNode* exonFeature,
                               int version) const;
    void setExonMappingVersion(const string& exonId,
                               vector<FeatureNode*> exonFeatures) const;
    void setExonsMappingVersions(ExonIdExonMap& exonIdExonMap) const;
    void setGeneMappingVersion(FeatureNode* geneTree) const;
    
    public:
    /* last minute fix-ups */
    void polishGene(FeatureNode* geneTree) const;
};

#endif
