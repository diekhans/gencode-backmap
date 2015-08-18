#include "featureMapper.hh"
#include "featureMapping.hh"
#include "gxf.hh"
#include "remapStatus.hh"
#include "pslMapping.hh"
#include "frame.hh"


/* 
 * create an feature for a full or partially mapped feature.
 */
void FeatureMapper::mkMappedFeature(const GxfFeature* feature,
                                    const PslCursor& srcPslCursor,
                                    const PslCursor& mappedPslCursor,
                                    int length,
                                    FeatureMapping* featureMapping) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int mappedTStart = mappedPslCursor.getTPosStrand('+');
    int mappedTEnd = mappedTStart + length;

    // add mapped feature. but don't update id or attributes now
    featureMapping->addMapped(
        gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                          mappedTStart, mappedTEnd, feature->fScore,
                          string(0, pslQStrand(srcPslCursor.getPsl())),
                          frame.toPhaseStr(), feature->fAttrs));
}

/* 
 * Create an feature for a full or partially unmapped feature.
 * The created feature is in source coordinates.
 * partIdx is used to make ID unique if split.
 */
void FeatureMapper::mkUnmappedFeature(const GxfFeature* feature,
                                      const PslCursor& srcPslCursor,
                                      const PslCursor& mappedPslCursor,
                                      int length,
                                      FeatureMapping* featureMapping) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int unmappedTStart = srcPslCursor.getTPosStrand('+');
    int unmappedTEnd = unmappedTStart + length;

    // add unmapped feature. but don't update id or attributes now
    featureMapping->addMapped(        
        gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                          unmappedTStart, unmappedTEnd, feature->fScore, feature->fStrand,
                          frame.toPhaseStr(), feature->fAttrs));
}

/*
 * Map one part of an feature.  Cursors are updated
 */
void FeatureMapper::mapFeaturePart(const GxfFeature* feature,
                                   PslCursor& srcPslCursor,
                                   PslCursor& mappedPslCursor,
                                   FeatureMapping* featureMapping) {
    assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
    if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
        // deleted region; length is minimum of different between starts in feature and how much is left in the feature
        int length = min(mappedPslCursor.getQPos()-srcPslCursor.getQPos(), srcPslCursor.getBlockLeft());
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length, featureMapping);
        srcPslCursor = srcPslCursor.advance(length);
    } else {
        // mapped region; length is the minimum left in either block
        int length = min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
        mkMappedFeature(feature, srcPslCursor, mappedPslCursor, length, featureMapping);
        srcPslCursor = srcPslCursor.advance(length);
        mappedPslCursor = mappedPslCursor.advance(length);
    }
}

/*
 * Map a feature.  Cursors are updated and srcPslCursor will point to the
 * next feature.
 */
void FeatureMapper::mapFeature(const GxfFeature* feature,
                               PslCursor& srcPslCursor,
                               PslCursor& mappedPslCursor,
                               FeatureMapping* featureMapping) {
    assert((feature->fEnd-feature->fStart)+1 == srcPslCursor.getBlockLeft());

    // note that source blocks can be merged in mapped block, but we don't merge
    // features.
    int srcPslFeatureQEnd = srcPslCursor.getQBlockEnd();
    while ((srcPslCursor.getQPos() < srcPslFeatureQEnd) && (not mappedPslCursor.atEnd())) {
        mapFeaturePart(feature, srcPslCursor, mappedPslCursor, featureMapping);
    }
    if (srcPslCursor.getQPos() < srcPslFeatureQEnd) {
        // unmapped at the end of feature; length is what is left over in this src block
        int length = srcPslCursor.getBlockLeft();
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length,  featureMapping);
        srcPslCursor = srcPslCursor.advance(length);
    }
    assert(srcPslCursor.getQPos() == srcPslFeatureQEnd);
}

/* mapped features */
void FeatureMapper::processMappedFeatures(const GxfFeatureVector& features,
                                          FeatureMappingSet* featureMappingSet) {
    PslCursor srcPslCursor(featureMappingSet->fPslMapping->fSrcPsl);
    PslCursor mappedPslCursor(featureMappingSet->fPslMapping->fMappedPsls[0]);
    for (int iFeature = 0; iFeature < features.size(); iFeature++) {
        FeatureMapping* featureMapping = featureMappingSet->add(features[iFeature]);
        mapFeature(features[iFeature], srcPslCursor, mappedPslCursor, featureMapping);
    }
}

/* process unmapped features, either sequence not in map, or no
 * mappings */
void FeatureMapper::processUnmappedFeatures(const GxfFeatureVector& features,
                                            FeatureMappingSet* featureMappingSet) {
    for (int iFeature = 0; iFeature < features.size(); iFeature++) {
        FeatureMapping* featureMapping = featureMappingSet->add(features[iFeature]);
        featureMapping->addUnmapped(features[iFeature]->clone());
    }
}

/* Map features of a transcript. Takes ownership of pslMapping object.
 * pslMapping will be NULL if source not in mapping alignments or when
 * indirect mappings can't be done because initial mapping is
 * deleted. Must explicitly pass srcSeqInMapping since we don't want to
 * flag indirect mapping features as not in mapping if the single features
 * are just not mapped */
FeatureMappingSet* FeatureMapper::mapFeatures(const GxfFeatureVector& features,
                                              PslMapping* pslMapping,
                                              bool srcSeqInMapping) {
    FeatureMappingSet* featureMappingSet = new FeatureMappingSet(pslMapping, srcSeqInMapping);
    if ((pslMapping == NULL) or not pslMapping->haveMappings()) {
        processUnmappedFeatures(features, featureMappingSet);
    } else {
        processMappedFeatures(features, featureMappingSet);
    }
    return featureMappingSet;
}

