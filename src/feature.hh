/* 
 * GENCODE features as a tree.
 */
#ifndef feature_hh
#define feature_hh
#include <assert.h>
#include "gxfRecord.hh"
#include "remapStatus.hh"



// FIXME: add the substituted target stuff caused this to step into hacky

/* Remap status attribute name */
extern const string REMAP_STATUS_ATTR;

/* Attribute name used for original id before remap */
extern const string REMAP_ORIGINAL_ID_ATTR;

/* Attribute name used for original location before remap */
extern const string REMAP_ORIGINAL_LOCATION_ATTR;

/* Attribute name for count of mappings, set on transcripts or genes */
extern const string REMAP_NUM_MAPPINGS_ATTR;

/* Attribute name for target of mapping */
extern const string REMAP_TARGET_STATUS_ATTR;

/* Attribute indicating target gene was substituted due to not being
 * able to map gene. Value is target version */
extern const string REMAP_SUBSTITUTED_MISSING_TARGET_ATTR;

class Feature;

/* Vector of Feature objects */
class FeatureVector: public vector<Feature*> {
public:
    /* sort the vector in a predictable order.  This is not necessary what
     * will be in the GxF file by GENCODE conventions. */
    void sort();
};

/**
 * Tree container for a GxfFeature object and children
 */
class Feature: public GxfFeature {
    private:
    Feature* fParent;
    FeatureVector fChildren;
    RemapStatus fRemapStatus;
    TargetStatus fTargetStatus;
    int fNumMappings;    // Number of location feature was mapped too.  Not set for all node types.

    private:
    int getTranscriptExonSize() const;
    int getOverlapAmount(const Feature* other) const;
    int countExonOverlap(const Feature* exon1,
                         const Feature* trans2) const;
    float getMaxTranscriptSimilarity(const Feature* gene2,
                                     const Feature* trans1,
                                     bool manualOnlyTranscripts) const;

    public:

    Feature(const string& seqid, const string& source, const string& type,
            int start, int end, const string& score, const string& strand,
            const string& phase, const AttrVals& attrs):
        GxfFeature(seqid, source, type, start, end, score, strand, phase, attrs),
        fParent(NULL),
        fRemapStatus(REMAP_STATUS_NONE),
        fTargetStatus(TARGET_STATUS_NA),
        fNumMappings(0) {
    }

    virtual ~Feature() {
        for (size_t i = 0; i < fChildren.size(); i++) {
            delete fChildren[i];
        }
    }

    /* clone the feature, not the tree */
    virtual Feature* clone() const {
        return new Feature(fSeqid, fSource, fType, fStart, fEnd, fScore, fStrand, fPhase, fAttrs);
    }
    
    /* clone tree */
    Feature* cloneTree() const;

    /* accessors */
    Feature* getParent() {
        return fParent;
    }
    FeatureVector& getChildren() {
        return fChildren;
    }
    const FeatureVector& getChildren() const {
        return fChildren;
    }
    Feature* getChild(int iChild) {
        return fChildren[iChild];
    }
    const Feature* getChild(int iChild) const {
        return fChildren[iChild];
    }
    RemapStatus getRemapStatus() const {
        return fRemapStatus;
    }
    TargetStatus getTargetStatus() const {
        return fTargetStatus;
    }

    int getNumMappings() const {
        return fNumMappings;
    }
    
    /* is this a gene? */
    bool isGene() const {
        return (fType == GxfFeature::GENE);
    }

    /* is this a transcript? */
    bool isTranscript() const {
        return (fType == GxfFeature::TRANSCRIPT);
    }

    /* is this an exon? */
    bool isExon() const {
        return (fType == GxfFeature::EXON);
    }

    /* is this a gene or transcript */
    bool isGeneOrTranscript() const {
        return isGene() or isTranscript();
    }

    /* get exon similarity */
    float getExonSimilarity(const Feature* trans2) const;

    /* get the maximum transcript similarity for a gene */
    float getMaxTranscriptSimilarity(const Feature* gene2,
                                     bool manualOnlyTranscripts=false) const;

    /* is ensembl small non-coding gene */
    bool isAutomaticSmallNonCodingGene() const;

    /* is this an automatic annotation? */
    bool isAutomatic() const {
        return fSource == GxfFeature::SOURCE_ENSEMBL;
    }

    /* is this an pseudogene annotation (excluding polymorphic)? */
    bool isPseudogene() const {
        const string& biotype = getTypeBiotype();
        return (biotype != "polymorphic_pseudogene")
            and (biotype.find("pseudogene") != biotype.npos);
    }

    bool anyChildWithRemapStatus(unsigned remapStatusSet) const;
    bool allChildWithRemapStatus(unsigned remapStatusSet) const;

    /* recursively get a list features matching the specified filter */
    void getMatching(FeatureVector& hits,
                     function<bool(const Feature*)>(filter)) const {
        if (filter(this)) {
            // FIXME: need const and non-const versions
            hits.push_back(const_cast<Feature*>(this));
        }
        for (int i = 0; i < fChildren.size(); i++) {
            fChildren[i]->getMatching(hits, filter);
        }
    }
    
    
    /* add a child node, linking up parents */
    void addChild(Feature* node) {
        assert(node->fParent == NULL);
        fChildren.push_back(node);
        node->fParent = this;
    }

    /* set remap status to specified value */
    void setRemapStatus(RemapStatus remapStatus) {
        fRemapStatus = remapStatus;
    }

    /* recursively set the remap status. */
    void rsetRemapStatus(RemapStatus remapStatus);

    /* Set the target status. Not recursive */
    void setTargetStatus(TargetStatus targetStatus);

    /* Recursively set the target status. */
    void rsetTargetStatus(TargetStatus targetStatus);

    /* recursively set the target status attribute */
    void rsetTargetStatusAttr();

    /* recursively set the target status attribute node. */
    void rsetSubstitutedMissingTargetAttr(const string& targetVersion);

    /* assign number of mappings */
    void setNumMappings(int numMappings) {
        fNumMappings = numMappings;
    }
    
    /* recursively print for debugging */
    void dump(ostream& fh) const;

    /* set the remap number of mappings attribute on this node  */
    void setNumMappingsAttr();

    /* recursively set the remap status attribute */
    void rsetRemapStatusAttr();

    /* depth-first output */
    void write(ostream& fh) const;
};

/* function to create a new feature based on type */
Feature* featureFactory(const string& seqid, const string& source, const string& type,
                        int start, int end, const string& score, const string& strand,
                        const string& phase, const AttrVals& attrs);

/* Get a base id, deleting the version, if it exists.
 */
static inline string getBaseId(const string& id) {
    size_t idot = id.find_last_of('.');
    return (idot == string::npos) ? id : id.substr(0, idot);
}

/* 
 * Get the id with mapping version (_N) removed, if it exists.
 */
static inline string getPreMappedId(const string& id) {
    size_t iun = id.find_last_of('_');
    if (iun == string::npos) {
        return id;
    } else {
        return id.substr(0, iun);
    }
}

/* 
 * Determine id an id has a mapping version (_N).
 */
static inline bool hasMappingVersion(const string& id) {
    size_t iun = id.find_last_of('_');
    return (iun != string::npos);
}

/* 
 * get the mapping version, or 0 if none
 */
static inline int getMappingVersion(const string& id) {
    size_t iun = id.find_last_of('_');
    if (iun == string::npos) {
        return 0;
    } else {
        return stringToInt(id.substr(iun+1));
    }
}

#endif
