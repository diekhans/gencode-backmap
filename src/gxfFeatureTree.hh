/*
 * Tree structure use to store genes.
 */
#ifndef gxfFeatureTree_hh
#define gxfFeatureTree_hh
#include <assert.h>
#include <functional>
#include "gxf.hh"

/**
 * Tree container for a GxfFeature object and children
 */
class GxfFeatureNode {
    public:
    const GxfFeature* fFeature;
    class GxfFeatureNode* fParent;
    vector<const GxfFeatureNode*> fChildren;

    GxfFeatureNode(const GxfFeature* feature):
        fFeature(feature),
        fParent(NULL) {
    }

    ~GxfFeatureNode() {
        delete fFeature;
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
    
    /* depth-first output */
    void write(ostream& fh) const {
        fGene->write(fh);
    }
};

#endif
