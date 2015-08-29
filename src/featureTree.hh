/*
 * Tree structure use to store genes.
 */
#ifndef featureTree_hh
#define featureTree_hh
#include <assert.h>
#include <functional>
#include "gxf.hh"
#include "remapStatus.hh"

/* Remap status attribute name */
extern const string REMAP_STATUS_ATTR;

/* Attribute name used for original id before remap */
extern const string REMAP_ORIGINAL_ID_ATTR;


/**
 * Tree container for a GxfFeature object and children
 */
class FeatureNode {
    public:
    GxfFeature* fFeature;
    FeatureNode* fParent;
    vector<FeatureNode*> fChildren;
    RemapStatus fRemapStatus;
    GxfFeatureVector fMappedFeatures;
    GxfFeatureVector fUnmappedFeatures;
    GxfFeatureVector fAllOutputFeatures;   // use for debugging, as it tracks order added.

    public:
    FeatureNode(GxfFeature* feature):
        fFeature(feature),
        fParent(NULL),
        fRemapStatus(REMAP_STATUS_NONE) {
    }

    ~FeatureNode() {
        delete fFeature;
        for (size_t i = 0; i < fMappedFeatures.size(); i++) {
            delete fMappedFeatures[i];
        }
        for (size_t i = 0; i < fUnmappedFeatures.size(); i++) {
            delete fUnmappedFeatures[i];
        }
        for (size_t i = 0; i < fChildren.size(); i++) {
            delete fChildren[i];
        }
    }

    /* recursively get a list features matching the specified filter */
    void getMatching(GxfFeatureVector& hits,
                     function<bool(const GxfFeature*)>(filter)) const {
        if (filter(fFeature)) {
            hits.push_back(fFeature);
        }
        for (int i = 0; i < fChildren.size(); i++) {
            fChildren[i]->getMatching(hits, filter);
        }
    }
    
    
    /* add a child node, linking up parents */
    void addChild(FeatureNode* node) {
        assert(node->fParent == NULL);
        fChildren.push_back(node);
        node->fParent = this;
    }

    /* add a mapped and take ownership */
    void addMapped(GxfFeature* mappedFeature) {
        fMappedFeatures.push_back(mappedFeature);
        fAllOutputFeatures.push_back(mappedFeature);
    }

    /* add a unmapped and take ownership */
    void addUnmapped(GxfFeature* unmappedFeature) {
        fUnmappedFeatures.push_back(unmappedFeature);
        fAllOutputFeatures.push_back(unmappedFeature);
    }

    /* compute the remap status of the feature. srcSeqInMapping
     * indicates of the srcSequence qs in the genomic map */
    RemapStatus calcRemapStatus(bool srcSeqInMapping) const;

    /* set the remap status, if remapStatus is REMAP_STATUS_NONE, calculate
     * it */
    void setRemapStatus(bool srcSeqInMapping,
                        RemapStatus remapStatus=REMAP_STATUS_NONE) {
        if (remapStatus == REMAP_STATUS_NONE) {
            fRemapStatus = calcRemapStatus(srcSeqInMapping);
        } else {
            fRemapStatus = remapStatus;
        }
    }
    
    /* print node for debugging */
    void dumpNode(ostream& fh) const;

    /* recursively print for debugging */
    void dump(ostream& fh) const;

    /* recursively set the remap status attribute */
    void setRemapStatusAttr();

    /* depth-first output */
    void write(ostream& fh) const;
};

/**
 * group a genes records together in a tree.
 */
class FeatureTree {
    private:
    void queueRecords(GxfParser *gxfParser,
                      GxfRecordVector& gxfRecords) const;
    FeatureNode* findGff3Parent(FeatureNode* geneTreeLeaf,
                                const GxfFeature* gxfFeature) const;
    FeatureNode* loadGff3GeneRecord(GxfFeature* gxfFeature,
                                    FeatureNode* geneTreeLeaf) const;
    const string& getGtfParentType(const string& featureType) const;
    FeatureNode* findGtfParent(FeatureNode* geneTreeLeaf,
                                  const GxfFeature* gxfFeature) const;
    FeatureNode* loadGtfGeneRecord(GxfFeature* gxfFeature,
                                   FeatureNode* geneTreeLeaf) const;    
    bool loadGeneRecord(GxfParser *gxfParser,
                        GxfRecord* gxfRecord,
                        FeatureNode* geneTreeRoot,
                        FeatureNode*& geneTreeLeaf,
                        GxfRecordVector& queuedRecords) const;
    FeatureNode* loadGene(GxfParser *gxfParser,
                          GxfFeature* geneFeature);
    public:
    FeatureNode* fGene;

    /* constructor */
    FeatureTree(GxfParser *gxfParser,
                GxfFeature* geneFeature);

    /* Destructor */
    ~FeatureTree();

    /* print for debugging */
    void dump(ostream& fh) const;

    /* depth-first output */
    void write(ostream& fh) const {
        fGene->write(fh);
    }
};

#endif
