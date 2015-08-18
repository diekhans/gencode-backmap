/*
 * Mapping of features, regardless of type
 */
#ifndef featureMapper_hh
#define featureMapper_hh

class GxfFeature;
class GxfFeatureVector;
class PslCursor;
class FeatureMapping;
class FeatureMappingSet;
class PslMapping;

/**
 * Static class to map a set of features of a transcript
 */
class FeatureMapper {
    private:
    static void mkMappedFeature(const GxfFeature* feature,
                                const PslCursor& srcPslCursor,
                                const PslCursor& mappedPslCursor,
                                int length,
                                FeatureMapping* featureMapping);
    static void mkUnmappedFeature(const GxfFeature* feature,
                                  const PslCursor& srcPslCursor,
                                  const PslCursor& mappedPslCursor,
                                  int length,
                                  FeatureMapping* featureMapping);
    static void mapFeaturePart(const GxfFeature* feature,
                               PslCursor& srcPslCursor,
                               PslCursor& mappedPslCursor,
                               FeatureMapping* featureMapping);
    static void mapFeature(const GxfFeature* feature,
                           PslCursor& srcPslCursor,
                           PslCursor& mappedPslCursor,
                           FeatureMapping* featureMapping);
    static void processMappedFeatures(const GxfFeatureVector& features,
                                      FeatureMappingSet* featureMappingSet);
    static void processUnmappedFeatures(const GxfFeatureVector& features,
                                        FeatureMappingSet* featureMappingSet);
    
    public:
    /* Map features of a transcript. Takes ownership of pslMapping object.
     * pslMapping will be NULL if source not in mapping alignments or when
     * indirect mappings can't be done because initial mapping is
     * deleted. Must explicitly pass srcSeqInMapping since we don't want to
     * flag indirect mapping features as not in mapping if the single features
     * are just not mapped */
    static FeatureMappingSet* mapFeatures(const GxfFeatureVector& features,
                                          PslMapping* pslMapping,
                                          bool srcSeqInMapping);
};
#endif
