/*
 * Tree structure use to store genes.
 */
#ifndef gxfFeatureTree_hh
#define gxfFeatureTree_hh
#include <assert.h>
#include <functional>
#include "gxf.hh"
#include "remapStatus.hh"

/**
 * Tree container for a GxfFeature object and children
 */
class GxfFeatureNode {
    public:
    const GxfFeature* fFeature;
    GxfFeatureNode* fParent;
    vector<GxfFeatureNode*> fChildren;
    RemapStatus fRemapStatus;
    GxfFeatureVector fMappedFeatures;
    GxfFeatureVector fUnmappedFeatures;

    private:
    GxfFeatureVector makeCombinedMappedUnmapped() const;

    public:
    GxfFeatureNode(const GxfFeature* feature):
        fFeature(feature),
        fParent(NULL),
        fRemapStatus(REMAP_STATUS_NONE) {
    }

    ~GxfFeatureNode() {
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
    void addChild(GxfFeatureNode* node) {
        assert(node->fParent == NULL);
        fChildren.push_back(node);
        node->fParent = this;
    }

    /* add a mapped and take ownership */
    void addMapped(const GxfFeature* mappedFeature) {
        fMappedFeatures.push_back(mappedFeature);
    }

    /* add a unmapped and take ownership */
    void addUnmapped(const GxfFeature* unmappedFeature) {
        fUnmappedFeatures.push_back(unmappedFeature);
    }

    /* compute the remap status of the feature. srcSeqInMapping
     * indicates of the srcSequence qs in the genomic map */
    RemapStatus calcRemapStatus(bool srcSeqInMapping) const;

    /* set the remap status */
    void setRemapStatus(bool srcSeqInMapping) {
        fRemapStatus = calcRemapStatus(srcSeqInMapping);
    }
    
    /* print node for debugging */
    void dumpNode(ostream& fh) const;

    /* recursively print for debugging */
    void dump(ostream& fh) const;

    /* depth-first output */
    void write(ostream& fh) const;
};

/**
 * group a genes records together in a tree.
 */
class GxfFeatureTree {
    private:
    void queueRecords(GxfParser *gxfParser,
                      GxfRecordVector& gxfRecords) const;
    GxfFeatureNode* findGff3Parent(GxfFeatureNode* geneTreeLeaf,
                                   const GxfFeature* gxfFeature) const;
    GxfFeatureNode* loadGff3GeneRecord(const GxfFeature* gxfFeature,
                                       GxfFeatureNode* geneTreeLeaf) const;
    const string& getGtfParentType(const string& featureType) const;
    GxfFeatureNode* findGtfParent(GxfFeatureNode* geneTreeLeaf,
                                  const GxfFeature* gxfFeature) const;
    GxfFeatureNode* loadGtfGeneRecord(const GxfFeature* gxfFeature,
                                      GxfFeatureNode* geneTreeLeaf) const;    
    bool loadGeneRecord(GxfParser *gxfParser,
                        const GxfRecord* gxfRecord,
                        GxfFeatureNode* geneTreeRoot,
                        GxfFeatureNode*& geneTreeLeaf,
                        GxfRecordVector& queuedRecords) const;
    GxfFeatureNode* loadGene(GxfParser *gxfParser,
                             const GxfFeature* geneFeature);
    public:
    GxfFeatureNode* fGene;

    /* constructor */
    GxfFeatureTree(GxfParser *gxfParser,
                   const GxfFeature* geneFeature);

    /* Destructor */
    ~GxfFeatureTree();

    /* print for debugging */
    void dump(ostream& fh) const;

    /* depth-first output */
    void write(ostream& fh) const {
        fGene->write(fh);
    }
};

#endif
