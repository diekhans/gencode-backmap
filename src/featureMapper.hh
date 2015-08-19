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
class GxfFeatureNode;
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
                                GxfFeatureNode* featureNode);
    static void mkUnmappedFeature(const GxfFeature* feature,
                                  const PslCursor& srcPslCursor,
                                  const PslCursor& mappedPslCursor,
                                  int length,
                                  GxfFeatureNode* featureNode);
    static void mapFeaturePart(const GxfFeature* feature,
                               PslCursor& srcPslCursor,
                               PslCursor& mappedPslCursor,
                               GxfFeatureNode* featureNode);
    static void mapFeature(const GxfFeature* feature,
                           PslCursor& srcPslCursor,
                           PslCursor& mappedPslCursor,
                           GxfFeatureNode* featureNode);
    static void processMappedFeature(GxfFeatureNode* featureNode,
                                     const PslMapping* pslMapping);
    static void processUnmappedFeature(GxfFeatureNode* featureNode,
                                       bool srcSeqInMapping);
    
    public:
    /* Map a single feature though an alignment of that feature.  The
     * pslMapping object will be NULL if source is not in mapping alignments
     * or when indirect mappings can't be done because initial mapping is
     * deleted.  Fill in mapped and unmapped arrays in featureNode. */
    static bool map(GxfFeatureNode* featureNode,
                    const PslMapping* pslMapping,
                    bool srcSeqInMapping);

    /* Map as single, bounding feature, like a gene or transcript record.
     * it's range is covered by contained ranges.  Omit new ranges if unmapped.
     * This puts all the adding of features in this one class.*/
    static void mapBounding(GxfFeatureNode* featureNode, bool srcSeqInMapping,
                            const string& targetSeqid="", int targetStart=-1, int targetEnd=-1, char targetStrand='.');
                       
};
#endif
