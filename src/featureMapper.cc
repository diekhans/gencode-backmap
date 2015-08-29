#include "featureMapper.hh"
#include "gxf.hh"
#include "featureTree.hh"
#include "remapStatus.hh"
#include "pslMapping.hh"
#include "frame.hh"
#include <iostream>

// FIXME: tmp
#define debug 0

/* 
 * create an feature for a full or partially mapped feature.
 */
void FeatureMapper::mkMappedFeature(const GxfFeature* feature,
                                    const PslCursor& srcPslCursor,
                                    const PslCursor& mappedPslCursor,
                                    int length,
                                    FeatureNode* featureNode) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int mappedTStart = mappedPslCursor.getTPosStrand('+');
    int mappedTEnd = mappedTStart + length;

    if (debug) {
        cerr << "mkMappedFeature: " << off << ": " << mappedTStart << "-" << mappedTEnd << endl;
    }

    // add mapped feature. but don't update id or attributes now
    featureNode->addMapped(
        gxfFeatureFactory(feature->getFormat(), string(mappedPslCursor.getPsl()->tName),
                          feature->fSource, feature->fType,
                          mappedTStart+1, mappedTEnd, feature->fScore,
                          charToString(pslQStrand(srcPslCursor.getPsl())),
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
                                      FeatureNode* featureNode) {
    int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
    Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

    // GxF genomic coordinates are always plus strand
    int unmappedTStart = srcPslCursor.getTPosStrand('+');
    int unmappedTEnd = unmappedTStart + length;

    if (debug) {
        cerr << "mkUnmappedFeature: " << off << ": " << unmappedTStart << "-" << unmappedTEnd << endl;
    }

    // add unmapped feature. but don't update id or attributes now
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
        if (debug) {
            cerr << "  mapFeaturePart: unmapped:  " << length << " " << srcPslCursor.toString() << " =>" << mappedPslCursor.toString() << endl;
        }
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length, featureNode);
        srcPslCursor = srcPslCursor.advance(length);
    } else {
        // mapped region; length is the minimum left in either block
        int length = min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
        if (debug) {
            cerr << "  mapFeaturePart: mapped:  " << length << " " << srcPslCursor.toString() << " =>" << mappedPslCursor.toString() << endl;
        }
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
    assert((feature->fEnd-feature->fStart)+1 == srcPslCursor.getBlockLeft());
    if (debug) {
        cerr << "mapFeature: " << feature->toString() << endl;
    }
    // note that source blocks can be merged in mapped block, but we don't merge
    // features.
    int srcPslFeatureQEnd = srcPslCursor.getQBlockEnd();
    while ((srcPslCursor.getQPos() < srcPslFeatureQEnd) && (not mappedPslCursor.atEnd())) {
        mapFeaturePart(feature, srcPslCursor, mappedPslCursor, featureNode);
    }
    if (srcPslCursor.getQPos() < srcPslFeatureQEnd) {
        // unmapped at the end of feature; length is what is left over in this src block
        int length = srcPslCursor.getBlockLeft();
        if (debug) {
            cerr << "  mapFeaturePart: unend:  " << length << " " << srcPslCursor.toString() << " =>" << mappedPslCursor.toString() << endl;
        }
        mkUnmappedFeature(feature, srcPslCursor, mappedPslCursor, length,  featureNode);
        srcPslCursor = srcPslCursor.advance(length);
    }
    assert(srcPslCursor.getQPos() == srcPslFeatureQEnd);
}

/* mapped features */
void FeatureMapper::processMappedFeature(FeatureNode* featureNode,
                                         const PslMapping* pslMapping) {
    PslCursor srcPslCursor(pslMapping->fSrcPsl);
    PslCursor mappedPslCursor(pslMapping->fMappedPsls[0]);
    mapFeature(featureNode->fFeature, srcPslCursor, mappedPslCursor, featureNode);
    featureNode->setRemapStatus(true);
}

/* process unmapped feature, either sequence not in map, or no
 * mappings */
void FeatureMapper::processUnmappedFeature(FeatureNode* featureNode,
                                           bool srcSeqInMapping) {
    featureNode->addUnmapped(featureNode->fFeature->clone());
    featureNode->setRemapStatus(srcSeqInMapping);
}

/* Map a single feature though an alignment of that feature.  The pslMapping
 * object will be NULL if source is not in mapping alignments or when indirect
 * mappings can't be done because initial mapping is deleted.  Fill in mapped
 * and unmapped arrays in featureNode. */
bool FeatureMapper::map(FeatureNode* featureNode,
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

/* Map as single, bounding feature, like a gene or transcript record.  it's
 * range is covered by contained ranges.  Omit new ranges if unmapped.
 */
void FeatureMapper::mapBounding(FeatureNode* featureNode, bool srcSeqInMapping,
                                const string& targetSeqid, int targetStart, int targetEnd, const string& targetStrand) {
    const GxfFeature* feature = featureNode->fFeature;
    if (targetStart >= 0) {
        featureNode->addMapped(
            gxfFeatureFactory(feature->getFormat(), targetSeqid,
                              feature->fSource, feature->fType,
                              targetStart+1, targetEnd, feature->fScore,
                              targetStrand, ".", feature->fAttrs));
    } else {
        featureNode->addUnmapped(featureNode->fFeature->clone());
    }
    featureNode->setRemapStatus(srcSeqInMapping);
}

/* Determine if an ID should be split into multiple unique ids. */
bool FeatureMapper::shouldSplitIds(const GxfFeatureVector& outputFeatures) {
    // FIXME: add check for discontinious ids.
    return (outputFeatures.size() > 1) && outputFeatures[0]->hasAttr("ID");
}


/* Assign a new id for a split attribute, save original id in remap_original_id
 */
void FeatureMapper::splitId(GxfFeature* outputFeature,
                            int partIdx) {
    const string& id = outputFeature->getAttrValue("ID");
    outputFeature->getAttrs().add(AttrVal(REMAP_ORIGINAL_ID_ATTR, id));
    outputFeature->getAttrs().update(AttrVal("ID", id+"_"+toString(partIdx)));
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
    const AttrVal* parentAttr = child->getAttrs().find("Parent");
    if (parentAttr != NULL) {
        child->getAttrs().update(AttrVal(parentAttr->fName, parent->getAttrValue("ID"),
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
void FeatureMapper::forceToUnmapped(FeatureNode* featureNode,
                                    RemapStatus remapStatus) {
    featureNode->fMappedFeatures.free();
    featureNode->fUnmappedFeatures.free();
    featureNode->addUnmapped(featureNode->fFeature->clone());
    featureNode->setRemapStatus(false, remapStatus);
}
