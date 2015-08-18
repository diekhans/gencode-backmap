/*
 * Mapping of features, regardless of type.
 */
#ifndef featureMapping_hh
#define featureMapping_hh
#include "gxf.hh"
#include "remapStatus.hh"

class PslMapping;

/**
 * Object to intermediately store mapping of a feature before they are added
 * to the GxF trees, since we need to look down the tree to set bounds and
 * status and up the tree to set parent ids.  Also, determine status requires
 * all the mappings of a feature, so it's easier to do this in two passes.
 */
class FeatureMapping {
    public:
    const GxfFeature* fSrcFeature;       // source feature
    const bool fSrcSeqInMapping;         // is the source sequence in the mapping?
    GxfFeatureVector fMappedFeatures;    // features that were mapped
    GxfFeatureVector fUnmappedFeatures;  // features that were not mapped
    
    public:
    /* constructor, has ownership of mapped/unmapped features */
    FeatureMapping(const GxfFeature* srcFeature,
                   bool srcSeqInMapping):
        fSrcFeature(srcFeature), fSrcSeqInMapping(srcSeqInMapping) {
    }

    /* destructor */
    ~FeatureMapping() {
        fMappedFeatures.free();
        fUnmappedFeatures.free();
    }

    /* add a mapped and take ownership */
    void addMapped(const GxfFeature* mappedFeature) {
        assert(fSrcSeqInMapping);
        fMappedFeatures.push_back(mappedFeature);
    }

    /* add a unmapped and take ownership */
    void addUnmapped(const GxfFeature* unmappedFeature) {
        fUnmappedFeatures.push_back(unmappedFeature);
    }

    /* compute the remap status of the feature */
    RemapStatus calcRemapStatus() const {
        if (not fSrcSeqInMapping) {
            // couldn't even try mapping, chrom not in map
            return REMAP_STATUS_NO_SEQ_MAP;
        } else if (fMappedFeatures.size() == 0) {
            assert(fUnmappedFeatures.size() > 0);
            // nothing mapped
            return REMAP_STATUS_DELETED;
        } else if (fUnmappedFeatures.size() == 0) {
            // full mapped
            if (fMappedFeatures.size() == 1) {
                return REMAP_STATUS_FULL_CONTIG;
            } else {
                return REMAP_STATUS_FULL_FRAGMENT;
            }
        } else {
            // partially mapped
            if (fMappedFeatures.size() == 1) {
                return REMAP_STATUS_PARTIAL_CONTIG;
            } else {
                return REMAP_STATUS_PARTIAL_FRAGMENT;
            }
        }
    }
};

/* Features that were mapped as a set from a single alignment.
 */
class FeatureMappingSet: public vector<FeatureMapping*> {
    public:
    const PslMapping* fPslMapping;
    const bool fSrcSeqInMapping;         // is the source sequence in the mapping?
    
    /* Constructor. Takes ownership of pslMapping object.  pslMapping will be
     * NULL if source not in mapping alignments or when indirect mappings
     * can't be done because initial mapping is deleted. Must explicitly pass
     * srcSeqInMapping since we don't want to flag indirect mapping features
     * as not in mapping if the single features are just not mapped */
    FeatureMappingSet(const PslMapping* pslMapping,
                      bool srcSeqInMapping):
        fPslMapping(pslMapping),
        fSrcSeqInMapping(srcSeqInMapping) {
    }

    /* destructor */
    ~FeatureMappingSet() {
        for (int i = 0; i < size(); i++) {
            delete (*this)[i];
        }
    }

    /* add a new feature mapping, doesn't take ownership of srcFeature */
    FeatureMapping* add(const GxfFeature* srcFeature) {
        push_back(new FeatureMapping(srcFeature, fSrcSeqInMapping));
        return (*this)[size()-1];
    }
};


#endif
