/*
 * Post mapping corrections to feature tree.
 */
#ifndef featureTreePolish_hh
#define featureTreePolish_hh
#include <assert.h>
#include <map>
#include "featureTree.hh"
class AnnotationSet;

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

    const AnnotationSet* fPreviousMappedAnotations; // maybe NULL
    
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
    bool isRemapped(const FeatureNode* featureNode) const;
    const FeatureNode* getPrevMappedFeature(const FeatureNode* newFeature) const;
    int getFeatureMappingVersion(const FeatureNode* prevFeature,
                                 bool featureSame) const;
    bool compareNodeAttrVals(const FeatureNode* prevNode,
                             const FeatureNode* newNode,
                             const string& attrName) const;
    bool compareNodeAttrs(const FeatureNode* prevNode,
                          const FeatureNode* newNode,
                          const StringVector& attrNames) const;
    bool compareMappedNodes(const FeatureNode* prevNode,
                            const FeatureNode* newNode,
                            const StringVector& attrNames) const;
    bool compareGeneNodes(const FeatureNode* prevNode,
                          const FeatureNode* newNode) const;
    bool compareTranscriptNodes(const FeatureNode* prevNode,
                                const FeatureNode* newNode) const;
    bool compareOtherNodes(const FeatureNode* prevNode,
                           const FeatureNode* newNode) const;
    bool compareMappedTranscripts(const FeatureNode* prevTranscript,
                                  const FeatureNode* newTranscript) const;
    bool compareMappedTranscriptsDescendants(const FeatureNode* prevParent,
                                             const FeatureNode* newParent) const;
    void setMappingVersionInId(FeatureNode* featureNode,
                               const AttrVal* attr,
                               int mappingVersion) const;
    void setMappingVersion(FeatureNode* featureNode,
                           const string& idAttrName,
                           const string& havanaIdAttrName,
                           int mappingVersion) const;
    void recursiveSetMappingVersion(FeatureNode* featureNode,
                                    const string& idAttrName,
                                    const string& havanaIdAttrName,
                                    int mappingVersion) const;
    void recordTranscriptMappedExons(FeatureNode* transcriptTree,
                                     ExonIdExonMap& exonIdExonMap) const;
    bool setTranscriptMappingVersion(FeatureNode* transcriptTree) const;
    void setExonMappingVersion(FeatureNode* exonFeature,
                               int mappingVersion) const;
    void setExonMappingVersion(const string& exonId,
                               vector<FeatureNode*> exonFeatures) const;
    void setExonsMappingVersions(ExonIdExonMap& exonIdExonMap) const;
    void setGeneMappingVersion(FeatureNode* geneTree) const;
    
    public:

    /* constructor */
    FeatureTreePolish(const AnnotationSet* previousMappedAnotations):
        fPreviousMappedAnotations(previousMappedAnotations) {
    }
    
    /* last minute fix-ups */
    void polishGene(FeatureNode* geneTree) const;
};

#endif
