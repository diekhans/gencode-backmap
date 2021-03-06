/*
 * Mapping of features, regardless of type.
 */
#ifndef featureMapper_hh
#define featureMapper_hh
#include <string>
#include "remapStatus.hh"
#include "featureTree.hh"
class TransMappedFeature;
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
    static FeatureNode* mkMappedFeature(const GxfFeature* feature,
                                        const PslCursor& srcPslCursor,
                                        const PslCursor& mappedPslCursor,
                                        int length);
    static FeatureNode* mkUnmappedFeature(const GxfFeature* feature,
                                          const PslCursor& srcPslCursor,
                                          const PslCursor& mappedPslCursor,
                                          int length);
    static void mapFeaturePart(const GxfFeature* feature,
                               PslCursor& srcPslCursor,
                               PslCursor& mappedPslCursor,
                               TransMappedFeature& transMappedFeature);
    static void mapFeature(const GxfFeature* feature,
                           PslCursor& srcPslCursor,
                           PslCursor& mappedPslCursor,
                           TransMappedFeature& transMappedFeature);
    static bool shouldSplitIds(const FeatureNodeVector& features);
    static void splitId(FeatureNode* feature,
                        int partIdx);
    static void splitIds(FeatureNodeVector& feature);
    static void splitIds(TransMappedFeature& transMappedFeature);
    static void processMappedFeature(const FeatureNode* feature,
                                     const PslMapping* pslMapping,
                                     TransMappedFeature& transMappedFeature);
    static void processUnmappedFeature(const FeatureNode* feature,
                                       TransMappedFeature& transMappedFeature);
    static FeatureNode* findContaining(FeatureNodeVector& parentFeatures,
                                       FeatureNode* childFeature);
    static void updateParent(FeatureNodeVector& parentFeatures,
                             FeatureNode* childFeature);
    static void updateParents(FeatureNodeVector& parentFeatures,
                              FeatureNodeVector& childFeatures);
    public:
    /* Map a single feature though an alignment of that feature.  The
     * pslMapping object will be NULL if source is not in mapping alignments
     * or when indirect mappings can't be done because initial mapping is
     * deleted. */
    static TransMappedFeature map(const FeatureNode* feature,
                               const PslMapping* pslMapping);

    /* update Parent id for mapped or unmapped, if needed. Link Feature
     * objects. */
    static void updateParent(FeatureNode* parentFeature,
                             FeatureNode* childFeature);

    /* validate parents and update Parent id for mapped and unmapped, in
     * needed. */
    static void updateParents(TransMappedFeature& parentFeatures,
                              TransMappedFeature& childFeatures);

    /* Map as single, bounding feature, like a gene or transcript record.
     * it's range is covered by contained ranges.  Omit new ranges if
     * unmapped.
     */
    static FeatureNode* mapBounding(const FeatureNode* feature,
                                    const string& targetSeqid="",
                                    int targetStart=-1,
                                    int targetEnd=-1,
                                    const string& targetStrand=".");
};
#endif
