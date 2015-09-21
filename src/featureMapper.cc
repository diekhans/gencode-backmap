#include "featureMapper.hh"
#include "gxf.hh"
#include "featureTree.hh"
#include "remapStatus.hh"
#include "pslMapping.hh"
#include "pslCursor.hh"
#include "frame.hh"
#include <iostream>

/* 
 * create an feature for a full or partially mapped feature.
 */
void FeatureMapper::mkMappedFeature(const GxfFeature* feature,
                                    const PslCursor& srcPslCursor,
                                    const PslCursor& mappedPslCursor,
                                    int length,
                                    FeatureNode* featureNode) {
    assert(length > 0);
    int off = srcPslCursor.getQPos() - srcPslCursor.getQStart();
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int mappedTStart, mappedTEnd;
    mappedPslCursor.getTRangeStrand('+', length, &mappedTStart, &mappedTEnd);

    // add mapped feature. but don't update id now
    GxfFeature* mappedFeature =
        gxfFeatureFactory(feature->getFormat(), string(mappedPslCursor.getPsl()->tName),
                          feature->fSource, feature->fType,
                          mappedTStart+1, mappedTEnd, feature->fScore,
                          charToString(pslTStrand(mappedPslCursor.getPsl())),
                          frame.toPhaseStr(), feature->fAttrs);
    featureNode->addMapped(mappedFeature);

    // save original coordinates for this region
    int srcTStart, srcTEnd;
    srcPslCursor.getTRangeStrand('+', length, &srcTStart, &srcTEnd);
    assert((feature->fStart-1 <= srcTStart) and (srcTEnd <= feature->fEnd));
    string origLocation = feature->fSeqid + ":" + feature->fStrand
        + ":" + toString(srcTStart+1) + "-" + toString(srcTEnd);
    mappedFeature->getAttrs().add(AttrVal(REMAP_ORIGINAL_LOCATION_ATTR, origLocation));
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
                                      FeatureNode* featureNode) {
    assert(length > 0);
    int off = srcPslCursor.getQPos() - srcPslCursor.getQStart();
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int unmappedTStart, unmappedTEnd;
    srcPslCursor.getTRangeStrand('+', length, &unmappedTStart, &unmappedTEnd);
    assert((feature->fStart-1 <= unmappedTStart) and (unmappedTEnd <= feature->fEnd));

    // add unmapped feature. but don't update id now. 
    featureNode->addUnmapped(
        gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                          unmappedTStart+1, unmappedTEnd, feature->fScore, feature->fStrand,
                          frame.toPhaseStr(), feature->fAttrs));
}

/*
 * Map one part of an feature.  Cursors are updated
 */
void FeatureMapper::mapFeaturePart(const GxfFeature* feature,
                                   PslCursor& srcPslCursor,
                                   PslCursor& mappedPslCursor,
                                   FeatureNode* featureNode) {
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
                               FeatureNode* featureNode) {
    assert(sameString(srcPslCursor.getPsl()->qName, mappedPslCursor.getPsl()->qName));
    assert(pslQStrand(srcPslCursor.getPsl()) == pslQStrand(mappedPslCursor.getPsl()));
    assert(srcPslCursor.getPsl()->qSize == mappedPslCursor.getPsl()->qSize);
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
void FeatureMapper::processMappedFeature(FeatureNode* featureNode,
                                         const PslMapping* pslMapping) {
    PslCursor srcPslCursor(pslMapping->fSrcPsl);
    PslCursor mappedPslCursor(pslMapping->fMappedPsl);
    mapFeature(featureNode->fFeature, srcPslCursor, mappedPslCursor, featureNode);
}

/* process unmapped feature, either sequence not in map, or no
 * mappings */
void FeatureMapper::processUnmappedFeature(FeatureNode* featureNode) {
    featureNode->addUnmapped(featureNode->fFeature->clone());
}

/* Map a single feature though an alignment of that feature.  The pslMapping
 * object will be NULL if source is not in mapping alignments or when indirect
 * mappings can't be done because initial mapping is deleted.  Fill in mapped
 * and unmapped arrays in featureNode. */
bool FeatureMapper::map(FeatureNode* featureNode,
                        const PslMapping* pslMapping) {
    if ((pslMapping == NULL) or not pslMapping->haveMappings()) {
        processUnmappedFeature(featureNode);
        return false;
    } else {
        processMappedFeature(featureNode, pslMapping);
        return true;
    }
}

/* compute status from children */
RemapStatus FeatureMapper::remapStatusFromChildren(FeatureNode* featureNode) {
    RemapStatus remapStatus = REMAP_STATUS_NONE;
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        remapStatus = remapStatusChildUpdate(remapStatus,
                                             featureNode->fChildren[i]->fRemapStatus);
    }
    return remapStatus;
}

/* Map as single, bounding feature, like a gene or transcript record.  it's
 * range is covered by contained ranges.  Omit new ranges if unmapped.
 * the total number of mappings.
 */
void FeatureMapper::mapBounding(FeatureNode* featureNode,
                                const string& targetSeqid,
                                int targetStart,
                                int targetEnd,
                                const string& targetStrand) {
    const GxfFeature* feature = featureNode->fFeature;
    if (targetStart >= 0) {
        GxfFeature* mappedFeature =
            gxfFeatureFactory(feature->getFormat(), targetSeqid,
                              feature->fSource, feature->fType,
                              targetStart+1, targetEnd, feature->fScore,
                              targetStrand, ".", feature->fAttrs);
        featureNode->addMapped(mappedFeature);
    } else {
        featureNode->addUnmapped(featureNode->fFeature->clone());
    }
}

/* Determine if an ID should be split into multiple unique ids. */
bool FeatureMapper::shouldSplitIds(const GxfFeatureVector& outputFeatures) {
    // FIXME: add check for discontinious ids.
    return (outputFeatures.size() > 1) && outputFeatures[0]->hasAttr(GxfFeature::ID_ATTR);
}


/* Assign a new id for a split attribute, save original id in remap_original_id
 */
void FeatureMapper::splitId(GxfFeature* outputFeature,
                            int partIdx) {
    const string& id = outputFeature->getAttrValue(GxfFeature::ID_ATTR);
    outputFeature->getAttrs().add(AttrVal(REMAP_ORIGINAL_ID_ATTR, id));
    outputFeature->getAttrs().update(AttrVal(GxfFeature::ID_ATTR, id+"_"+toString(partIdx)));
}


/* If a feature was split into multiple mapped or unmapped parts, assign new ids.
 * Add remap_original_id attribute.
 * FIXME: this should not split discontinuous features.
 */
void FeatureMapper::splitIds(GxfFeatureVector& outputFeatures) {
    for (int iPart = 0; iPart < outputFeatures.size(); iPart++) {
        splitId(outputFeatures[iPart], iPart);
    }
}

/* If a feature was split into multiple mapped and unmapped parts, assign new ids.
 * Add remap_original_id attribute.
 */
void FeatureMapper::splitIds(FeatureNode* featureNode) {
    if (shouldSplitIds(featureNode->fMappedFeatures)) {
        splitIds(featureNode->fMappedFeatures);
    }
    if (shouldSplitIds(featureNode->fUnmappedFeatures)) {
        splitIds(featureNode->fUnmappedFeatures);
    }
}

/* containing parent feature in a list, or error if not found */
GxfFeature* FeatureMapper::findContaining(GxfFeature* child,
                                          GxfFeatureVector& parentParts) {
    for (int i = 0; i < parentParts.size(); i++) {
        if ((parentParts[i]->fStart <= child->fStart)
            and (child->fEnd <= parentParts[i]->fEnd)) {
            return parentParts[i];
        }
    }
    throw invalid_argument("BUG: containing parent not found: " + child->toString());
}

/* validate parents and update Parent id for mapped or unmapped, if needed. */
void FeatureMapper::updateParent(GxfFeature* child,
                                 GxfFeatureVector& parentParts) {
    GxfFeature* parent = findContaining(child, parentParts);  // validates even if not used
    const AttrVal* parentAttr = child->getAttrs().find(GxfFeature::PARENT_ATTR);
    if (parentAttr != NULL) {
        child->getAttrs().update(AttrVal(parentAttr->fName, parent->getAttrValue(GxfFeature::ID_ATTR),
                                         parentAttr->fQuoted));
    }
}

/* validate parents and update Parent id for mapped or unmapped, in needed. */
void FeatureMapper::updateParents(GxfFeatureVector& childParts,
                                  GxfFeatureVector& parentParts) {
    for (int i = 0; i < childParts.size(); i++) {
        updateParent(childParts[i], parentParts);
    }
}

/* validate parents of mapped and unppaed, Parent id if needed. */
void FeatureMapper::updateParents(FeatureNode* featureNode,
                                  FeatureNode* parentNode) {
    updateParents(featureNode->fMappedFeatures,
                  parentNode->fMappedFeatures);
    updateParents(featureNode->fUnmappedFeatures,
                  parentNode->fUnmappedFeatures);
}

/*
 * Recursively update node ids and parent links if this is a GFF3.  Also
 * verifies that child node is contained withing a parent.  This is all that
 * is done for GTF.  Parent is passed down rather than found by pointer in
 * structure so we don't attempt to update above the tree we started with.
 */
void FeatureMapper::updateIds(FeatureNode* featureNode,
                              FeatureNode* parentNode) {
    splitIds(featureNode);
    if (parentNode != NULL) {
        updateParents(featureNode, parentNode);
    }
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        updateIds(featureNode->fChildren[i], featureNode);
    }
}

/* Convert a feature to full unmapped, removing any partial unmapped features.
 * Used when transcript conflicts withing a gene are found.
 */
void FeatureMapper::forceToUnmapped(FeatureNode* featureNode) {
    featureNode->fMappedFeatures.free();
    featureNode->fUnmappedFeatures.free();
    featureNode->addUnmapped(featureNode->fFeature->clone());
}
