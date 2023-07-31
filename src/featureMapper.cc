#include "featureMapper.hh"
#include "featureTree.hh"
#include "remapStatus.hh"
#include "pslMapping.hh"
#include "pslCursor.hh"
#include "frame.hh"
#include <iostream>

/* 
 * create an feature for a full or partially mapped feature.
 */
FeatureNode* FeatureMapper::mkMappedFeature(const GxfFeature* feature,
                                            const PslCursor& srcPslCursor,
                                            const PslCursor& mappedPslCursor,
                                            int length) {
    assert(length > 0);
    int off = srcPslCursor.getQPos() - srcPslCursor.getQStart();
    Frame frame(Frame::fromPhaseStr(feature->getPhase()).incr(off));

    // GxF genomic coordinates are always plus strand
    int mappedTStart, mappedTEnd;
    mappedPslCursor.getTRangeStrand('+', length, &mappedTStart, &mappedTEnd);

    // add mapped feature. but don't update id now
    FeatureNode* mappedFeature
        = FeatureNode::factory(string(mappedPslCursor.getPsl()->tName),
                               feature->getSource(), feature->getType(),
                               mappedTStart+1, mappedTEnd, feature->getScore(),
                               charToString(pslTStrand(mappedPslCursor.getPsl())),
                               frame.toPhaseStr(), feature->getAttrs());
                          
    // save original coordinates for this region
    int srcTStart, srcTEnd;
    srcPslCursor.getTRangeStrand('+', length, &srcTStart, &srcTEnd);
    assert((feature->getStart()-1 <= srcTStart) and (srcTEnd <= feature->getEnd()));
    string origLocation = feature->getSeqid() + ":" + feature->getStrand()
        + ":" + toString(srcTStart+1) + "-" + toString(srcTEnd);
    mappedFeature->getAttrs().add(AttrVal(REMAP_ORIGINAL_LOCATION_ATTR, origLocation));
    return mappedFeature;
}

/* 
 * Create an feature for a full or partially unmapped feature.
 * The created feature is in source coordinates.
 * partIdx is used to make ID unique if split.
 */
FeatureNode* FeatureMapper::mkUnmappedFeature(const GxfFeature* feature,
                                              const PslCursor& srcPslCursor,
                                              const PslCursor& mappedPslCursor,
                                              int length) {
    assert(length > 0);
    int off = srcPslCursor.getQPos() - srcPslCursor.getQStart();
    Frame frame(Frame::fromPhaseStr(feature->getPhase()).incr(off));

    // GxF genomic coordinates are always plus strand
    int unmappedTStart, unmappedTEnd;
    srcPslCursor.getTRangeStrand('+', length, &unmappedTStart, &unmappedTEnd);
    assert((feature->getStart()-1 <= unmappedTStart) and (unmappedTEnd <= feature->getEnd()));

    // add unmapped feature. but don't update id now. 
    FeatureNode* unmappedFeature =
        FeatureNode::factory(feature->getSeqid(), feature->getSource(), feature->getType(),
                             unmappedTStart+1, unmappedTEnd, feature->getScore(), feature->getStrand(),
                             frame.toPhaseStr(), feature->getAttrs());
    return unmappedFeature;
}

/*
 * Map one part of an feature.  Cursors are updated
 */
void FeatureMapper::mapFeaturePart(const GxfFeature* feature,
                                   PslCursor& srcPslCursor,
                                   PslCursor& mappedPslCursor,
                                   TransMappedFeature& transMappedFeature) {
    assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
    if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
        // deleted region; length is minimum of different between starts in feature and how much is left in the feature
        int length = min(mappedPslCursor.getQPos()-srcPslCursor.getQPos(), srcPslCursor.getBlockLeft());
        transMappedFeature.addUnmapped(mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length));
        srcPslCursor = srcPslCursor.advance(length);
    } else {
        // mapped region; length is the minimum left in either block
        int length = min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
        transMappedFeature.addMapped(mkMappedFeature(feature, srcPslCursor, mappedPslCursor, length));
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
                               TransMappedFeature& transMappedFeature) {
    assert(sameString(srcPslCursor.getPsl()->qName, mappedPslCursor.getPsl()->qName));
    assert(pslQStrand(srcPslCursor.getPsl()) == pslQStrand(mappedPslCursor.getPsl()));
    assert(srcPslCursor.getPsl()->qSize == mappedPslCursor.getPsl()->qSize);
    assert((feature->getEnd() - feature->getStart())+1 == srcPslCursor.getBlockLeft());

    // note that source blocks can be merged in mapped block, but we don't merge
    // features.
    int srcPslFeatureQEnd = srcPslCursor.getQBlockEnd();
    while ((srcPslCursor.getQPos() < srcPslFeatureQEnd) && (not mappedPslCursor.atEnd())) {
        mapFeaturePart(feature, srcPslCursor, mappedPslCursor, transMappedFeature);
    }
    if (srcPslCursor.getQPos() < srcPslFeatureQEnd) {
        // unmapped at the end of feature; length is what is left over in this src block
        int length = srcPslCursor.getBlockLeft();
        transMappedFeature.addUnmapped(mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length));
        srcPslCursor = srcPslCursor.advance(length);
    }
    assert(srcPslCursor.getQPos() == srcPslFeatureQEnd);
}

/* Determine if an ID should be split into multiple unique ids. */
bool FeatureMapper::shouldSplitIds(const FeatureNodeVector& features) {
    // FIXME: add check for discontinious ids.
    return (features.size() > 1) && features[0]->hasAttr(GxfFeature::ID_ATTR);
}

/* Assign a new id for a split attribute, save original id in remap_original_id
 */
void FeatureMapper::splitId(FeatureNode* feature,
                            int partIdx) {
    const string& id = feature->getAttrValue(GxfFeature::ID_ATTR);
    feature->getAttrs().add(AttrVal(REMAP_ORIGINAL_ID_ATTR, id));
    feature->getAttrs().update(AttrVal(GxfFeature::ID_ATTR, id + "_" + toString(partIdx)));
}

/* If a feature was split into multiple mapped or unmapped parts, assign new ids.
 * Add remap_original_id attribute.
 * FIXME: this should not split discontinuous features.
 */
void FeatureMapper::splitIds(FeatureNodeVector& features) {
    for (int iPart = 0; iPart < features.size(); iPart++) {
        splitId(features[iPart], iPart);
    }
}

/* split ids if needed */
void FeatureMapper::splitIds(TransMappedFeature& transMappedFeature) {
    if (shouldSplitIds(transMappedFeature.mapped)) {
        splitIds(transMappedFeature.mapped);
    }
    if (shouldSplitIds(transMappedFeature.unmapped)) {
        splitIds(transMappedFeature.unmapped);
    }
}

/* mapped features */
void FeatureMapper::processMappedFeature(const FeatureNode* feature,
                                         const PslMapping* pslMapping,
                                         TransMappedFeature& transMappedFeature) {
    PslCursor srcPslCursor(pslMapping->getSrcPsl());
    PslCursor mappedPslCursor(pslMapping->getBestMappedPsl());
    mapFeature(feature->getGxfFeature(), srcPslCursor, mappedPslCursor, transMappedFeature);
}

/* process unmapped feature, either sequence not in map, or no
 * mappings */
void FeatureMapper::processUnmappedFeature(const FeatureNode* feature,
                                           TransMappedFeature& transMappedFeature) {
    transMappedFeature.addUnmapped(feature->cloneFeature());
}

/* Map a single feature though an alignment of that feature.  The pslMapping
 * object will be NULL if source is not in mapping alignments or when indirect
 * mappings can't be done because initial mapping is deleted. */
TransMappedFeature FeatureMapper::map(const FeatureNode* feature,
                                      const PslMapping* pslMapping) {
    TransMappedFeature transMappedFeature(feature);
    if ((pslMapping == NULL) or not pslMapping->haveMappings()) {
        processUnmappedFeature(feature, transMappedFeature);
    } else {
        processMappedFeature(feature, pslMapping, transMappedFeature);
    }
    splitIds(transMappedFeature);
    return transMappedFeature;
}

/* containing parent feature in a list, or error if not found */
FeatureNode* FeatureMapper::findContaining(FeatureNodeVector& parentFeatures,
                                           FeatureNode* childFeature) {
    for (int i = 0; i < parentFeatures.size(); i++) {
        if ((parentFeatures[i]->getStart() <= childFeature->getStart())
            and (childFeature->getEnd() <= parentFeatures[i]->getEnd())) {
            return parentFeatures[i];
        }
    }
    throw invalid_argument("BUG: containing parent not found: " + childFeature->toString());
}

/* update Parent id for mapped or unmapped, if needed. Link Feature
 * objects. */
void FeatureMapper::updateParent(FeatureNode* parentFeature,
                                 FeatureNode* childFeature) {
    parentFeature->getChildren().push_back(childFeature);
    const AttrVal* parentAttr = childFeature->getAttrs().find(GxfFeature::PARENT_ATTR);
    if (parentAttr != NULL) {
        childFeature->getAttrs().update(AttrVal(parentAttr->getName(),
                                                parentFeature->getAttrValue(GxfFeature::ID_ATTR)));
    }
}

/* validate parents and update Parent id for mapped or unmapped, if
 * needed. Link Feature objects. */
void FeatureMapper::updateParent(FeatureNodeVector& parentFeatures,
                                 FeatureNode* childFeature) {
    FeatureNode* parentFeature = findContaining(parentFeatures, childFeature);
    updateParent(parentFeature, childFeature);
}

/* validate parents and update Parent id for mapped or unmapped, in needed. */
void FeatureMapper::updateParents(FeatureNodeVector& parentFeatures,
                                  FeatureNodeVector& childFeatures) {
    for (int i = 0; i < childFeatures.size(); i++) {
        updateParent(parentFeatures, childFeatures[i]);
    }
}

/* validate parents and update Parent id for mapped and unmapped, in needed. */
void FeatureMapper::updateParents(TransMappedFeature& parentFeatures,
                                  TransMappedFeature& childFeatures) {
    if (parentFeatures.mapped.size() > 0) {
        updateParents(parentFeatures.mapped, childFeatures.mapped);
    }
    if (parentFeatures.unmapped.size() > 0) {
        updateParents(parentFeatures.unmapped, childFeatures.unmapped);
    }
}

/* Map as single, bounding feature, like a gene or transcript record.  it's
 * range is covered by contained ranges.  Omit new ranges if unmapped.
 * the total number of mappings.
 */
FeatureNode* FeatureMapper::mapBounding(const FeatureNode* feature,
                                        const string& targetSeqid,
                                        int targetStart,
                                        int targetEnd,
                                        const string& targetStrand) {
    // FIXME: could parent update be here?
    if (targetStart >= 0) {
        return FeatureNode::factory(targetSeqid, feature->getSource(), feature->getType(),
                                    targetStart+1, targetEnd, feature->getScore(),
                                    targetStrand, ".", feature->getAttrs());
    } else {
        // clone only feature, not tree.
        return feature->cloneFeature();
    }
}
