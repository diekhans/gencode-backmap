/**
 * Structure for storing results during processing
 */
#ifndef resultFeatures_hh
#define resultFeatures_hh
#include "feature.hh"
#include "remapStatus.hh"
#include <iostream>

/**
 * Trees resulting from a map.
 */
class ResultFeatures {
    private:
    /* do any child belond to the specified status */
    bool anyChildWithRemapStatus(unsigned remapStatusSet) const {
        return ((mapped != NULL) && mapped->anyChildWithRemapStatus(remapStatusSet))
            || ((unmapped != NULL) && unmapped->anyChildWithRemapStatus(remapStatusSet));
    }

    /* do all child have belong to the specified status set */
    bool allChildWithRemapStatus(unsigned remapStatusSet) const {
        return ((mapped == NULL) || ((mapped != NULL) && mapped->allChildWithRemapStatus(remapStatusSet)))
            && ((unmapped == NULL) || ((unmapped != NULL) && unmapped->allChildWithRemapStatus(remapStatusSet)));
    }


    public:
    const Feature* src;  /// not owned
    Feature* mapped;
    Feature* unmapped;
    Feature* target;    // substituted from target

    /* constructor */
    ResultFeatures(const Feature* src = NULL,
                   Feature* mapped = NULL,
                   Feature* unmapped = NULL):
        src(src), mapped(mapped), unmapped(unmapped), target(NULL) {
    }

    /* free data all */
    void free() {
        delete mapped;
        delete unmapped;
        delete target;
    }
    /* free mapped */
    void freeMapped() {
        delete mapped;
        mapped = NULL;
    }

    /* free unmapped */
    void freeUnmapped() {
        delete unmapped;
        unmapped = NULL;
    }

    /* get remap status from either mappedm unmapped, or target. */
    RemapStatus getRemapStatus() const {
        if (mapped != NULL) {
            return mapped->getRemapStatus();
        } else if (unmapped != NULL) {
            return unmapped->getRemapStatus();
        } else if (target != NULL) {
            return target->getRemapStatus();
        } else {
            return REMAP_STATUS_DELETED;
        }
    }

    /* recursively set the remap status on mapped and unmapped */
    void rsetRemapStatus(RemapStatus remapStatus) {
        if (mapped != NULL) {
            mapped->rsetRemapStatus(remapStatus);
        }
        if (unmapped != NULL) {
            unmapped->rsetRemapStatus(remapStatus);
        }
    }

    /* get target status from either mapped, unmapped, or target */
    TargetStatus getTargetStatus() const {
        if (mapped != NULL) {
            return mapped->getTargetStatus();
        } else if (unmapped != NULL) {
            return unmapped->getTargetStatus();
        } else if (target != NULL) {
            return target->getTargetStatus();
        } else {
            return TARGET_STATUS_LOST;
        }
    }

    /* set target status on mapped and unmapped trees */
    void setTargetStatus(TargetStatus targetStatus) {
        if (mapped != NULL) {
            mapped->setTargetStatus(targetStatus);
        }
        if (unmapped != NULL) {
            unmapped->setTargetStatus(targetStatus);
        }
    }

    /* recursively set target status on trees */
    void rsetTargetStatus(TargetStatus targetStatus) {
        if (mapped != NULL) {
            mapped->rsetTargetStatus(targetStatus);
        }
        if (unmapped != NULL) {
            unmapped->rsetTargetStatus(targetStatus);
        }
    }

    /* recursively set the target status attribute */
    void rsetTargetStatusAttr() {
        if (mapped != NULL) {
            mapped->rsetTargetStatusAttr();
        }
        if (unmapped != NULL) {
            unmapped->rsetTargetStatusAttr();
        }
    }

    /* get number of mappings from either mapped or unmapped. */
    int getNumMappings() const {
        if (mapped != NULL) {
            return mapped->getNumMappings();
        } else if (unmapped != NULL) {
            return unmapped->getNumMappings();
        } else {
            return 0;
        }
    }

    /* set the remap number of mappings attribute on this node (not
     * recursive) */
    void setNumMappingsAttr() {
        if (mapped != NULL) {
            mapped->setNumMappingsAttr();
        }
        if (unmapped != NULL) {
            unmapped->setNumMappingsAttr();
        }
    }

    /* recursively set the remap status attribute */
    void rsetRemapStatusAttr() {
        if (mapped != NULL) {
            mapped->rsetRemapStatusAttr();
        }
        if (unmapped != NULL) {
            unmapped->rsetRemapStatusAttr();
        }
    }

    /* determine gene or transcript remap status from children.  This doesn't
     * handle GENE_CONFLICT or_GENE_SIZE_CHANGE, which are forced.
     */
    RemapStatus calcBoundingFeatureRemapStatus(bool srcSeqInMapping) const {
        assert(src->isGeneOrTranscript());
        if (not srcSeqInMapping) {
            return REMAP_STATUS_NO_SEQ_MAP;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_FULL_CONTIG)) {
            return REMAP_STATUS_FULL_CONTIG;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_DELETED)) {
            return REMAP_STATUS_DELETED;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_FULL_CONTIG|REMAP_STATUS_FULL_FRAGMENT)) {
            return REMAP_STATUS_FULL_FRAGMENT;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_FULL_CONTIG|REMAP_STATUS_FULL_FRAGMENT|REMAP_STATUS_PARTIAL|REMAP_STATUS_DELETED)) {
            return REMAP_STATUS_PARTIAL;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_GENE_CONFLICT)) {
            return REMAP_STATUS_GENE_CONFLICT;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_GENE_SIZE_CHANGE)) {
            return REMAP_STATUS_GENE_SIZE_CHANGE;
        }
        if (allChildWithRemapStatus(REMAP_STATUS_AUTOMATIC_GENE)) {
            return REMAP_STATUS_AUTOMATIC_GENE;
        }
        dump(cerr);
        throw logic_error("gene RemapStatus logic error");
    }

    /* set remap status on a bounding */
    void setBoundingFeatureRemapStatus(bool srcSeqInMapping) {
        RemapStatus remapStatus = calcBoundingFeatureRemapStatus(srcSeqInMapping);
        if (mapped != NULL) {
            mapped->setRemapStatus(remapStatus);
        }
        if (unmapped != NULL) {
            unmapped->setRemapStatus(remapStatus);
        }
    }

    /* print for debugging */
    void dump(ostream& fh) const {
        if (mapped != NULL) {
            fh << "@@@ mapped" << endl;
            mapped->dump(fh);
        }
        if (unmapped != NULL) {
            fh << "@@@ unmapped" << endl;
            unmapped->dump(fh);
        }
    }
};

/* vector of mapped resulting feature */
class ResultFeaturesVector: public vector<ResultFeatures> {
    public:
    bool haveMapped() const {
        for (int i = 0; i < size(); i++) {
            if (((*this)[i]).mapped != NULL) {
                return true;
            }
        }
        return false;
    }
    bool haveUnmapped() const {
        for (int i = 0; i < size(); i++) {
            if (((*this)[i]).unmapped != NULL) {
                return true;
            }
        }
        return false;
    }
};
/**
 * Set of mapped and unmapped features resulting from transmap.  A feature maybe split when mapped,
 * hence vectors.
 */
class TransMappedFeature {
    public:
    const Feature* src;
    FeatureVector mapped;
    FeatureVector unmapped;

    /* constructors */
    TransMappedFeature(const Feature* src = NULL):
        src(src) {
    }

    TransMappedFeature(ResultFeatures& resultFeatures) {
        src = resultFeatures.src;
        if (resultFeatures.mapped != NULL) {
            mapped.push_back(resultFeatures.mapped);
        }
        if (resultFeatures.unmapped != NULL) {
            unmapped.push_back(resultFeatures.unmapped);
        }
    }

    /** add a mapped node */
    void addMapped(Feature* feature) {
        mapped.push_back(feature);
    }

    /** add a unmapped node */
    void addUnmapped(Feature* feature) {
        unmapped.push_back(feature);
    }

    /* compute status for a transmapped feature.  this only looks at a single
     * level of mappings, not a tree. srcSeqInMapping indicates of the srcSequence
     * was in the genomic map */
    RemapStatus calcRemapStatus(bool srcSeqInMapping) const {
        if (not srcSeqInMapping) {
            // couldn't even try mapping, chrom not in map
            return REMAP_STATUS_NO_SEQ_MAP;
        } else if (mapped.size() == 0) {
            assert(unmapped.size() > 0);
            // nothing mapped
            return REMAP_STATUS_DELETED;
        } else if (unmapped.size() == 0) {
            // full mapped
            if (mapped.size() == 1) {
                return REMAP_STATUS_FULL_CONTIG;
            } else {
                return REMAP_STATUS_FULL_FRAGMENT;
            }
        } else {
            // partially mapped
            assert(mapped.size() > 0);
            assert(unmapped.size() > 0);
            return REMAP_STATUS_PARTIAL;
        }
    }

    /**
     * calculate the remap status and set it in the feature nodes
     */
    void setRemapStatus(bool srcSeqInMapping) {
        RemapStatus remapStatus = calcRemapStatus(srcSeqInMapping);
        for (int i = 0; i < mapped.size(); i++) {
            mapped[i]->setRemapStatus(remapStatus);
        }
        for (int i = 0; i < unmapped.size(); i++) {
            unmapped[i]->setRemapStatus(remapStatus);
        }
    }
};

#endif
