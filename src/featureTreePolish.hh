/*
 * Post mapping corrections to feature tree.
 */
#ifndef featureTreePolish_hh
#define featureTreePolish_hh
#include <assert.h>
#include <map>
#include "feature.hh"
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
    typedef map<int, FeatureVector> ExonNumExonMap;
    typedef ExonNumExonMap::iterator ExonNumExonMapIter;
    typedef ExonNumExonMap::const_iterator ExonNumExonMapConstIter;

    // use to add mapping version number to exon ids
    typedef map<string, FeatureVector> ExonIdExonMap;
    typedef ExonIdExonMap::iterator ExonIdExonMapIter;
    typedef ExonIdExonMap::const_iterator ExonIdExonMapConstIter;

    const AnnotationSet* fPreviousMappedAnotations; // maybe NULL
    
    void renumberExon(Feature* exon,
                      int exonNum,
                      ExonNumExonMap& exonNumExonMap) const;
    void renumberExons(Feature* transcript,
                       ExonNumExonMap& exonNumExonMap) const;
    Feature* findNewExon(Feature* feature,
                             int oldExonNum,
                             ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeature(Feature* feature,
                              ExonNumExonMap& exonNumExonMap) const;
    void renumberOtherFeatures(Feature* feature,
                               ExonNumExonMap& exonNumExonMap) const;
    void renumberTranscriptExons(Feature* transcript) const;
    void renumberGeneExons(Feature* gene) const;
    bool isRemapped(const Feature* feature) const;
    const Feature* getPrevMappedFeature(const Feature* newFeature) const;
    int getFeatureMappingVersion(const Feature* prevFeature,
                                 bool featureSame) const;
    bool compareAttrVal(const AttrVal* prevAttr,
                        const AttrVal* newAttr,
                        int iValue,
                        bool isIdAttr) const;
    bool compareAttrVals(const Feature* prevFeature,
                         const Feature* newFeature,
                         const string& attrName,
                         bool isIdAttr) const;
    bool compareAttrs(const Feature* prevFeature,
                      const Feature* newFeature,
                      const StringVector& attrNames,
                      bool isIdAttrs) const;
    bool compareAttrs(const Feature* prevFeature,
                      const Feature* newFeature,
                      const StringVector& attrNames,
                      const StringVector& idAttrNames) const;
    bool compareMappedFeatures(const Feature* prevFeature,
                               const Feature* newFeature,
                               const StringVector& attrNames,
                               const StringVector& idAttrNames) const;
    bool compareGeneFeatures(const Feature* prevFeature,
                             const Feature* newFeature) const;
    bool compareTranscriptFeatures(const Feature* prevFeature,
                                   const Feature* newFeature) const;
    bool compareOtherFeatures(const Feature* prevFeature,
                              const Feature* newFeature) const;
    bool compareMappedTranscripts(const Feature* prevTranscript,
                                  const Feature* newTranscript) const;
    bool compareMappedTranscriptsDescendants(const Feature* prevParent,
                                             const Feature* newParent) const;
    void setMappingVersionInId(Feature* feature,
                               const AttrVal* attr,
                               int mappingVersion) const;
    void setMappingVersion(Feature* feature,
                           const string& idAttrName,
                           const string& havanaIdAttrName,
                           int mappingVersion) const;
    void recursiveSetMappingVersion(Feature* feature,
                                    const string& idAttrName,
                                    const string& havanaIdAttrName,
                                    int mappingVersion) const;
    bool setTranscriptMappingVersion(Feature* transcript) const;
    void collectExons(const Feature* root,
                      ExonIdExonMap& exonIdExonMap) const;
    int getExonMappingVersion(FeatureVector& exonFeatures) const;
    void setExonMappingVersion(Feature* exonFeature,
                               int mappingVersion) const;
    void setExonMappingVersion(const string& exonId,
                               FeatureVector& exonFeatures,
                               FeatureVector* prevExonFeatures) const;
    void setExonsMappingVersions(const Feature* prevGene,
                                 Feature* gene) const;
    bool setTranscriptsMappingVersions(Feature* gene) const;
    void setGeneMappingVersion(Feature* gene) const;
    
    public:

    /* constructor */
    FeatureTreePolish(const AnnotationSet* previousMappedAnotations):
        fPreviousMappedAnotations(previousMappedAnotations) {
    }
    
    /* last minute fix-ups */
    void polishGene(Feature* gene) const;
};

#endif
