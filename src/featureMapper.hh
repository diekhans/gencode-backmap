/*
 * Mapping of features, regardless of type.
 */
#ifndef featureMapper_hh
#define featureMapper_hh
#include <string>
using namespace std;


class GxfFeature;
class GxfFeatureVector;
class PslCursor;
class FeatureNode;
class PslMapping;

/**
 * Class of static functions to map a set of features of a transcript, using
 * an mapping of the feature coordinates.
 */
class FeatureMapper {
    private:
    static void mkMappedFeature(const GxfFeature* feature,
                                const PslCursor& srcPslCursor,
                                const PslCursor& mappedPslCursor,
                                int length,
                                FeatureNode* featureNode);
    static void mkUnmappedFeature(const GxfFeature* feature,
                                  const PslCursor& srcPslCursor,
                                  const PslCursor& mappedPslCursor,
                                  int length,
                                  FeatureNode* featureNode);
    static void mapFeaturePart(const GxfFeature* feature,
                               PslCursor& srcPslCursor,
                               PslCursor& mappedPslCursor,
                               FeatureNode* featureNode);
    static void mapFeature(const GxfFeature* feature,
                           PslCursor& srcPslCursor,
                           PslCursor& mappedPslCursor,
                           FeatureNode* featureNode);
    static void processMappedFeature(FeatureNode* featureNode,
                                     const PslMapping* pslMapping);
    static void processUnmappedFeature(FeatureNode* featureNode,
                                       bool srcSeqInMapping);
    static bool shouldSplitIds(const GxfFeatureVector& outputFeatures);
    static void splitId(GxfFeature* outputFeature,
                        int partIdx);
    static void splitId(GxfFeatureVector& outputFeatures,
                        int partIdx);
    static void splitIds(GxfFeatureVector& outputFeatures);
    static void splitIds(FeatureNode* featureNode);
    static GxfFeature* findContaining(GxfFeature* child,
                                      GxfFeatureVector& parentParts);
    static void updateParent(GxfFeature* child,
                             GxfFeatureVector& parentParts);
    static void updateParents(GxfFeatureVector& childParts,
                              GxfFeatureVector& parentParts);
    static void updateParents(FeatureNode* featureNode,
                              FeatureNode* parentNode);

    public:
    /* Map a single feature though an alignment of that feature.  The
     * pslMapping object will be NULL if source is not in mapping alignments
     * or when indirect mappings can't be done because initial mapping is
     * deleted.  Fill in mapped and unmapped arrays in featureNode. */
    static bool map(FeatureNode* featureNode,
                    const PslMapping* pslMapping,
                    bool srcSeqInMapping);

    /* Map as single, bounding feature, like a gene or transcript record.
     * it's range is covered by contained ranges.  Omit new ranges if
     * unmapped.
     */
    static void mapBounding(FeatureNode* featureNode, bool srcSeqInMapping,
                            const string& targetSeqid="", int targetStart=-1, int targetEnd=-1, const string& targetStrand=".");

    /*
     * Recursively update node ids and parent links if this is a GFF3.  Also
     * verifies that child node is contained withing a parent.  This is all that
     * is done for GTF.  Parent is passed down rather than found by pointer in
     * structure so we don't attempt to update above the tree we started with.
     */
    static void updateIds(FeatureNode* featureNode,
                          FeatureNode* parentNode=NULL);
};
#endif
