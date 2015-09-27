/*
 * Tree structure use to store genes
 */
#include "featureTree.hh"
#include <ostream>
#include <algorithm> 

/* Remap status attribute name */
const string REMAP_STATUS_ATTR = "remap_status";

/* Attribute name used for original id before remap */
const string REMAP_ORIGINAL_ID_ATTR = "remap_original_id";

/* Attribute name used for original localization before remap */
const string REMAP_ORIGINAL_LOCATION_ATTR = "remap_original_location";

/* Attribute name for count of mappings, set on transcripts or genes */
const string REMAP_NUM_MAPPINGS_ATTR = "remap_num_mappings";

/* Attribute name for target of mapping */
const string REMAP_TARGET_STATUS_ATTR = "remap_target_status";

/* Attribute indicating target gene was substituted due to   */
const string REMAP_SUBSTITUTED_MISSING_TARGET_ATTR = "remap_substituted_missing_target";


/* set the remap number of mappings attribute on this node  */
void FeatureNode::setNumMappingsAttr() {
    AttrVal numMappingsAttr(REMAP_NUM_MAPPINGS_ATTR, toString(fNumMappings));
    for (int i = 0; i < fMappedFeatures.size(); i++) {
        fMappedFeatures[i]->getAttrs().add(numMappingsAttr);
    }
    for (int i = 0; i < fUnmappedFeatures.size(); i++) {
        fUnmappedFeatures[i]->getAttrs().add(numMappingsAttr);
    }
}

/* recursively set the remap status attribute */
void FeatureNode::setRemapStatusAttr() {
    AttrVal remapStatusAttr(REMAP_STATUS_ATTR, remapStatusToStr(fRemapStatus));
    for (int i = 0; i < fMappedFeatures.size(); i++) {
        fMappedFeatures[i]->getAttrs().add(remapStatusAttr);
    }
    for (int i = 0; i < fUnmappedFeatures.size(); i++) {
        fUnmappedFeatures[i]->getAttrs().add(remapStatusAttr);
    }
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->setRemapStatusAttr();
    }
}

/* recursively set the target status attribute */
void FeatureNode::setTargetStatusAttr() {
    AttrVal targetStatusAttr(REMAP_TARGET_STATUS_ATTR, targetStatusToStr(fTargetStatus));
    for (int i = 0; i < fMappedFeatures.size(); i++) {
        fMappedFeatures[i]->getAttrs().add(targetStatusAttr);
    }
    for (int i = 0; i < fUnmappedFeatures.size(); i++) {
        fUnmappedFeatures[i]->getAttrs().add(targetStatusAttr);
    }
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->setTargetStatusAttr();
    }
}

/* depth-first output */
void FeatureNode::write(ostream& fh) const {
    fh << fFeature->toString() << endl;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->write(fh);
    }
}


/* Return a list records, moving from gxfRecords vector  */
void GeneTree::queueRecords(GxfParser *gxfParser,
                            GxfRecordVector& gxfRecords) {
    for (size_t i = 0; i < gxfRecords.size(); i++) {
        gxfParser->push(gxfRecords[i]);
    }
    gxfRecords.clear();
}

/*
 * Find the parent for GFF3.
 */
FeatureNode* GeneTree::findGff3Parent(FeatureNode* geneTreeLeaf,
                                      const GxfFeature* gxfFeature) {
    const string& parentId = gxfFeature->getAttrValue(GxfFeature::PARENT_ATTR);
    FeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->getAttrValue(GxfFeature::ID_ATTR) != parentId)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node " +  parentId + " for " + gxfFeature->getAttrValue(GxfFeature::ID_ATTR) + " not found");
    }
    return parent;
}

/*
 * Process a FF3 record for a gene, which uses the explicit tree.
 * Return the new leaf node.
 */
FeatureNode* GeneTree::loadGff3GeneRecord(GxfFeature* gxfFeature,
                                          FeatureNode* geneTreeLeaf) {
    FeatureNode* parent = findGff3Parent(geneTreeLeaf, gxfFeature);
    FeatureNode* child = new FeatureNode(gxfFeature);
    parent->addChild(child);
    return child;
}

/* Get the desired type the parent of GTF feature
 * WARNING: this assumed that the hierarchy is
 * gene->transcript->{everything else}
 * FIXME: this is what GFF3 does, which might not be right.
 */
const string& GeneTree::getGtfParentType(const string& featureType) {
    assert(featureType != GxfFeature::GENE);
    if (featureType == GxfFeature::TRANSCRIPT) {
        return GxfFeature::GENE;
    } else {
        return GxfFeature::TRANSCRIPT;
    }
}

/*
 * Find the parent for a GTF record.  This is painful guess based on the
 * GENCODE file order and know how GENCODE is structures.
 */
FeatureNode* GeneTree::findGtfParent(FeatureNode* geneTreeLeaf,
                                     const GxfFeature* gxfFeature) {
    const string& parentType = getGtfParentType(gxfFeature->fType);
    FeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->fType != parentType)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node of type " + parentType + "  not found for type " + parent->fFeature->fType);
    }
    return parent;
}

/* compute the remap status of the feature. srcSeqInMapping
 * indicates of the srcSequence qs in the genomic map */
RemapStatus FeatureNode::calcRemapStatus(bool srcSeqInMapping) const {
    if (not srcSeqInMapping) {
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
        return REMAP_STATUS_PARTIAL;
    }
}

/* recursively determine the remap status */
void FeatureNode::setRemapStatusFromChildren(RemapStatus baseStatus) {
    fRemapStatus = baseStatus;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fRemapStatus = remapStatusChildUpdate(fRemapStatus, fChildren[i]->fRemapStatus);
    }
    
}

/* recursively determine the remap status */
void FeatureNode::recursiveCalcRemapStatus(bool srcSeqInMapping) {
    // compute leaves first
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->recursiveCalcRemapStatus(srcSeqInMapping);
    }

    // calculate for node and considering children
    setRemapStatusFromChildren(calcRemapStatus(srcSeqInMapping));
}

/* print node for debugging */
void FeatureNode::dumpNode(ostream& fh) const {
    const string status = remapStatusToStr(fRemapStatus);
    fh << "src" << "\t" << ((fFeature == NULL) ? "NULL" : fFeature->toString()) << endl;
    for (int i = 0; i < fAllOutputFeatures.size(); i++) {
        fh << (fMappedFeatures.contains(fAllOutputFeatures[i]) ? "mapped" : "unmapped") << "\t" << status << "\t" << fAllOutputFeatures[i]->toString() << endl;
    }
}

/* recursively print for debugging */
void FeatureNode::dump(ostream& fh) const {
    dumpNode(fh);
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->dump(fh);
    }
}

/*
 * Process a GTF record for a gene, which uses knowledge of
 * the GENCODE structure to reproduce the hierarchy.
 * Return the new leaf node.
 */
FeatureNode* GeneTree::loadGtfGeneRecord(GxfFeature* gxfFeature,
                                         FeatureNode* geneTreeLeaf) {
    FeatureNode* parent = findGtfParent(geneTreeLeaf, gxfFeature);
    FeatureNode* child = new FeatureNode(gxfFeature);
    parent->addChild(child);
    return child;
}

/*
 * Process a GxfRecord for a gene, return False if no more for this gene.
 */
bool GeneTree::loadGeneRecord(GxfParser *gxfParser,
                              GxfRecord* gxfRecord,
                              FeatureNode* geneTreeRoot,
                              FeatureNode*& geneTreeLeaf,
                              GxfRecordVector& queuedRecords) {
    if (instanceOf(gxfRecord, GxfLine)) {
        queuedRecords.push_back(gxfRecord);
        return true;
    } else {
        GxfFeature* gxfFeature = dynamic_cast<GxfFeature*>(gxfRecord);
        if (gxfFeature->fType == GxfFeature::GENE) {
            queuedRecords.push_back(gxfRecord); // next gene
            return false;
        } else {
            if (gxfParser->getGxfFormat() == GFF3_FORMAT) {
                geneTreeLeaf = loadGff3GeneRecord(gxfFeature, geneTreeLeaf);
            } else {
                geneTreeLeaf = loadGtfGeneRecord(gxfFeature, geneTreeLeaf);
            }
            return true;
        }
    }
}

/*
 * Load all records associated with a given gene.  Return non-feature and the
 * next gene to the queue to process.  This will causes comments in the middle
 * of genes to be moved to the end, but GENCODE doesn't put comments in genes.
 * This whole thing is annoying due to the lack of explicit structure in GTF.
 */
FeatureNode* GeneTree::loadGene(GxfParser *gxfParser,
                                GxfFeature* geneFeature) {
    assert(geneFeature->fType == GxfFeature::GENE);

    FeatureNode* geneTreeRoot = new FeatureNode(geneFeature);
    FeatureNode* geneTreeLeaf = geneTreeRoot;  // were we are currently working
    GxfRecordVector queuedRecords;
    GxfRecord* gxfRecord = NULL;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (not loadGeneRecord(gxfParser, gxfRecord, geneTreeRoot, geneTreeLeaf, queuedRecords)) {
            break;
        }
    }
    queueRecords(gxfParser, queuedRecords);
    return geneTreeRoot;
}

/* factory */
FeatureNode* GeneTree::geneTreeFactory(GxfParser *gxfParser,
                                       GxfFeature* geneFeature) {
    return loadGene(gxfParser, geneFeature);
}
