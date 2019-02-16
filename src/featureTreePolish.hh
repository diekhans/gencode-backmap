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
    typedef map<int, FeatureNodeVector> ExonNumExonMap;
    typedef ExonNumExonMap::iterator ExonNumExonMapIter;
    typedef ExonNumExonMap::const_iterator ExonNumExonMapConstIter;

    // use to add mapping version number to exon ids
    typedef map<string, FeatureNodeVector> ExonIdExonMap;
    typedef ExonIdExonMap::iterator ExonIdExonMapIter;
    typedef ExonIdExonMap::const_iterator ExonIdExonMapConstIter;

    const AnnotationSet* fPreviousMappedAnotations; // maybe NULL
    
    void renumberExon(FeatureNode* exon,
                      int exonNum,
                      ExonNumExonMap& exonNumExonMap) const;
    void renumberExons(FeatureNode* transcript,
                       ExonNumExonMap& exonNumExonMap) const;
    FeatureNode* findNewExon(FeatureNode* feature,
                             int oldExonNum,
                             ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeature(FeatureNode* feature,
                              ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeatures(FeatureNode* feature,
                               ExonNumExonMap& exonNumExonMap) const;
    void renumberTranscriptExons(FeatureNode* transcript) const;
    void renumberGeneExons(FeatureNode* gene) const;
    bool isRemapped(const FeatureNode* feature) const;
    const FeatureNode* getPrevMappedFeature(const FeatureNode* newFeature) const;
    int getFeatureMappingVersion(const FeatureNode* prevFeature,
                                 bool featureSame) const;
    bool compareAttrVal(const AttrVal* prevAttr,
                        const AttrVal* newAttr,
                        int iValue,
                        bool isIdAttr) const;
    bool compareAttrVals(const FeatureNode* prevFeature,
                         const FeatureNode* newFeature,
                         const string& attrName,
                         bool isIdAttr) const;
    bool compareAttrs(const FeatureNode* prevFeature,
                      const FeatureNode* newFeature,
                      const StringVector& attrNames,
                      bool isIdAttrs) const;
    bool compareAttrs(const FeatureNode* prevFeature,
                      const FeatureNode* newFeature,
                      const StringVector& attrNames,
                      const StringVector& idAttrNames) const;
    bool compareMappedFeatures(const FeatureNode* prevFeature,
                               const FeatureNode* newFeature,
                               const StringVector& attrNames,
                               const StringVector& idAttrNames) const;
    bool compareGeneFeatures(const FeatureNode* prevFeature,
                             const FeatureNode* newFeature) const;
    bool compareTranscriptFeatures(const FeatureNode* prevFeature,
                                   const FeatureNode* newFeature) const;
    bool compareOtherFeatures(const FeatureNode* prevFeature,
                              const FeatureNode* newFeature) const;
    bool compareMappedTranscripts(const FeatureNode* prevTranscript,
                                  const FeatureNode* newTranscript) const;
    bool compareMappedTranscriptsDescendants(const FeatureNode* prevParent,
                                             const FeatureNode* newParent) const;
    void setMappingVersionInId(FeatureNode* feature,
                               const AttrVal* attr,
                               int mappingVersion) const;
    void setMappingVersion(FeatureNode* feature,
                           const string& idAttrName,
                           const string& havanaIdAttrName,
                           int mappingVersion) const;
    void recursiveSetMappingVersion(FeatureNode* feature,
                                    const string& idAttrName,
                                    const string& havanaIdAttrName,
                                    int mappingVersion) const;
    bool setTranscriptMappingVersion(FeatureNode* transcript) const;
    void collectExons(const FeatureNode* root,
                      ExonIdExonMap& exonIdExonMap) const;
    int getExonMappingVersion(FeatureNodeVector& exonFeatures) const;
    void setExonMappingVersion(FeatureNode* exonFeature,
                               int mappingVersion) const;
    void setExonMappingVersion(const string& exonId,
                               FeatureNodeVector& exonFeatures,
                               FeatureNodeVector* prevExonFeatures) const;
    void setExonsMappingVersions(const FeatureNode* prevGene,
                                 FeatureNode* gene) const;
    bool setTranscriptsMappingVersions(FeatureNode* gene) const;
    void setGeneMappingVersion(FeatureNode* gene) const;
    
    public:

    /* constructor */
    FeatureTreePolish(const AnnotationSet* previousMappedAnotations):
        fPreviousMappedAnotations(previousMappedAnotations) {
    }
    
    /* last minute fix-ups */
    void polishGene(FeatureNode* gene) const;
};

#endif
