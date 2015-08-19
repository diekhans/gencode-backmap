#include "featureMapper.hh"
#include "gxf.hh"
#include "gxfFeatureTree.hh"
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
                                    GxfFeatureNode* featureNode) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int mappedTStart = mappedPslCursor.getTPosStrand('+');
    int mappedTEnd = mappedTStart + length;

    // add mapped feature. but don't update id or attributes now
    featureNode->addMapped(
        gxfFeatureFactory(feature->getFormat(), string(mappedPslCursor.getPsl()->tName),
                          feature->fSource, feature->fType,
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
                                      GxfFeatureNode* featureNode) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int unmappedTStart = srcPslCursor.getTPosStrand('+');
    int unmappedTEnd = unmappedTStart + length;

    // add unmapped feature. but don't update id or attributes now
    featureNode->addMapped(        
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
                                   GxfFeatureNode* featureNode) {
    assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
    if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
        // deleted region; length is minimum of different between starts in feature and how much is left in the feature
        int length = min(mappedPslCursor.getQPos()-srcPslCursor.getQPos(), srcPslCursor.getBlockLeft());
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length, featureNode);
        srcPslCursor = srcPslCursor.advance(length);
    } else {
        // mapped region; length is the minimum left in either block
        int length = min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
        mkMappedFeature(feature, srcPslCursor, mappedPslCursor, length, featureNode);
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
                               GxfFeatureNode* featureNode) {
    assert((feature->fEnd-feature->fStart)+1 == srcPslCursor.getBlockLeft());

    // note that source blocks can be merged in mapped block, but we don't merge
    // features.
    int srcPslFeatureQEnd = srcPslCursor.getQBlockEnd();
    while ((srcPslCursor.getQPos() < srcPslFeatureQEnd) && (not mappedPslCursor.atEnd())) {
        mapFeaturePart(feature, srcPslCursor, mappedPslCursor, featureNode);
    }
    if (srcPslCursor.getQPos() < srcPslFeatureQEnd) {
        // unmapped at the end of feature; length is what is left over in this src block
        int length = srcPslCursor.getBlockLeft();
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length,  featureNode);
        srcPslCursor = srcPslCursor.advance(length);
    }
    assert(srcPslCursor.getQPos() == srcPslFeatureQEnd);
}

/* mapped features */
void FeatureMapper::processMappedFeature(GxfFeatureNode* featureNode,
                                         const PslMapping* pslMapping) {
    PslCursor srcPslCursor(pslMapping->fSrcPsl);
    PslCursor mappedPslCursor(pslMapping->fMappedPsls[0]);
    mapFeature(featureNode->fFeature, srcPslCursor, mappedPslCursor, featureNode);
    featureNode->setRemapStatus(true);
}

/* process unmapped feature, either sequence not in map, or no
 * mappings */
void FeatureMapper::processUnmappedFeature(GxfFeatureNode* featureNode,
                                           bool srcSeqInMapping) {
    featureNode->addUnmapped(featureNode->fFeature->clone());
    featureNode->setRemapStatus(srcSeqInMapping);
}

/* Map a single feature though an alignment of that feature.  The pslMapping
 * object will be NULL if source is not in mapping alignments or when indirect
 * mappings can't be done because initial mapping is deleted.  Fill in mapped
 * and unmapped arrays in featureNode. */
bool FeatureMapper::map(GxfFeatureNode* featureNode,
                        const PslMapping* pslMapping,
                        bool srcSeqInMapping) {
    if ((pslMapping == NULL) or not pslMapping->haveMappings()) {
        processUnmappedFeature(featureNode, srcSeqInMapping);
        return false;
    } else {
        processMappedFeature(featureNode, pslMapping);
        return true;
    }
}

/* Map as single, bounding feature, like a gene or transcript record.
 * it's range is covered by contained ranges.  Omit new ranges if unmapped.
 * This puts all the adding of features in this one class.*/
void FeatureMapper::mapBounding(GxfFeatureNode* featureNode, bool srcSeqInMapping,
                                const string& targetSeqid, int targetStart, int targetEnd, char targetStrand) {
    const GxfFeature* feature = featureNode->fFeature;
    if (targetStart >= 0) {
        featureNode->addMapped(
            gxfFeatureFactory(feature->getFormat(), targetSeqid,
                              feature->fSource, feature->fType,
                              targetStart+1, targetEnd, feature->fScore,
                              string(0, targetStrand), ".", feature->fAttrs));
    } else {
        featureNode->addUnmapped(featureNode->fFeature->clone());
    }
    featureNode->setRemapStatus(srcSeqInMapping);
}
                       

