/*
 * Tree structure use to store genes.
 */
#ifndef gxfFeatureTree_hh
#define gxfFeatureTree_hh
#include <assert.h>
#include "gxf.hh"

/* status of remap of feature */
typedef enum {
    REMAP_STATUS_NONE,               // not set
    REMAP_STATUS_FULL_CONTIG,        // full remap, contiguous
    REMAP_STATUS_FULL_FRAGMENT,      // full remap, fragmented
    REMAP_STATUS_PARTIAL_CONTIG,     // partial remap, contiguous
    REMAP_STATUS_PARTIAL_FRAGMENT,   // partial remap, fragmented
    REMAP_STATUS_DELETED,            // completely deleted
    REMAP_STATUS_NO_SEQ_MAP          // no mapping alignments for sequence
} RemapStatus;

/* Attribute name used for remap status */
extern const string REMAP_ATTR_NAME;

/* convert a remap status to a AttrVal  */
const AttrVal& remapStatusToAttrVal(RemapStatus remapStatus);
   

/**
 * Tree container for a GxfFeature object and children
 */
class GxfFeatureNode {
    public:
    const GxfFeature* fFeature;
    const RemapStatus fRemapStatus;
    class GxfFeatureNode* fParent;
    vector<const GxfFeatureNode*> fChildren;

    GxfFeatureNode(const GxfFeature* feature,
                   RemapStatus remapStatus):
        fFeature(feature),
        fRemapStatus(remapStatus),
        fParent(NULL) {
    }

    ~GxfFeatureNode() {
        delete fFeature;
        for (size_t i = 0; i < fChildren.size(); i++) {
            delete fChildren[i];
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
