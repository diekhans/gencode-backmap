/*
 * Tree structure use to store genes
 */
#ifndef gxfFeatureTree_hh
#define gxfFeatureTree_hh
#include <assert.h>

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

    /* add a child node, linking up parents */
    void addChild(GxfFeatureNode* node) {
        assert(node->fParent == NULL);
        fChildren.push_back(node);
        node->fParent = this;
    }
};

#endif
