/*
 * Mapping of features, regardless of type.
 */
#ifndef featureMapper_hh
#define featureMapper_hh
#include <string>
#include "remapStatus.hh"
#include "feature.hh"
class TransMappedFeature;
using namespace std;


class GxfFeature;
class GxfFeatureVector;
class PslCursor;
class Feature;
class PslMapping;


/**
 * Class of static functions to map a set of features of a transcript, using
 * an mapping of the feature coordinates.
 */
class FeatureMapper {
    private:
    static Feature* mkMappedFeature(const GxfFeature* feature,
                                        const PslCursor& srcPslCursor,
                                        const PslCursor& mappedPslCursor,
                                        int length);
    static Feature* mkUnmappedFeature(const GxfFeature* feature,
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
    static bool shouldSplitIds(const FeatureVector& features);
    static void splitId(Feature* feature,
                        int partIdx);
    static void splitIds(FeatureVector& feature);
    static void splitIds(TransMappedFeature& transMappedFeature);
    static void processMappedFeature(const Feature* feature,
                                     const PslMapping* pslMapping,
                                     TransMappedFeature& transMappedFeature);
    static void processUnmappedFeature(const Feature* feature,
                                       TransMappedFeature& transMappedFeature);
    static Feature* findContaining(FeatureVector& parentFeatures,
                                       Feature* childFeature);
    static void updateParent(FeatureVector& parentFeatures,
                             Feature* childFeature);
    static void updateParents(FeatureVector& parentFeatures,
                              FeatureVector& childFeatures);
    public:
    /* Map a single feature though an alignment of that feature.  The
     * pslMapping object will be NULL if source is not in mapping alignments
     * or when indirect mappings can't be done because initial mapping is
     * deleted. */
    static TransMappedFeature map(const Feature* feature,
                               const PslMapping* pslMapping);

    /* update Parent id for mapped or unmapped, if needed. Link Feature
     * objects. */
    static void updateParent(Feature* parentFeature,
                             Feature* childFeature);

    /* validate parents and update Parent id for mapped and unmapped, in
     * needed. */
    static void updateParents(TransMappedFeature& parentFeatures,
                              TransMappedFeature& childFeatures);

    /* Map as single, bounding feature, like a gene or transcript record.
     * it's range is covered by contained ranges.  Omit new ranges if
     * unmapped.
     */
    static Feature* mapBounding(const Feature* feature,
                                const string& targetSeqid="",
                                int targetStart=-1,
                                int targetEnd=-1,
                                const string& targetStrand=".");
};
#endif
