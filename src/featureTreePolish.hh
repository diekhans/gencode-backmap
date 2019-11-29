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
    public:
    // use to map old exon numbers to new exons in a transcript
    typedef map<int, FeatureNodeVector> ExonNumExonMap;
    typedef ExonNumExonMap::iterator ExonNumExonMapIter;
    typedef ExonNumExonMap::const_iterator ExonNumExonMapConstIter;

    // use to add mapping version number to exon ids
    typedef map<string, FeatureNodeVector> ExonIdExonMap;
    typedef ExonIdExonMap::iterator ExonIdExonMapIter;
    typedef ExonIdExonMap::const_iterator ExonIdExonMapConstIter;

    private:
    const AnnotationSet* fPreviousMappedAnotations; // maybe NULL
    
    const FeatureNode* getPrevMappedFeature(const FeatureNode* newFeature) const;
    bool setTranscriptMappingVersion(FeatureNode* transcript) const;
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
