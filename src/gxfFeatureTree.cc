/*
 * Tree structure use to store genes
 */
#include "gxfFeatureTree.hh"
#include <ostream>
#include <algorithm> 

/* depth-first output */
void GxfFeatureNode::write(ostream& fh) const {
    fh << fFeature->toString() << endl;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->write(fh);
    }
}


/* Return a list records, moving from gxfRecords vector  */
void GxfFeatureTree::queueRecords(GxfParser *gxfParser,
                                  GxfRecordVector& gxfRecords) const {
    for (size_t i = 0; i < gxfRecords.size(); i++) {
        gxfParser->push(gxfRecords[i]);
    }
    gxfRecords.clear();
}

/*
 * Find the parent for GFF3.
 */
GxfFeatureNode* GxfFeatureTree::findGff3Parent(GxfFeatureNode* geneTreeLeaf,
                                               const GxfFeature* gxfFeature) const {
    const string& parentId = gxfFeature->getAttrValue("Parent");
    GxfFeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->getAttrValue("ID") != parentId)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node " +  parentId + " for " + gxfFeature->getAttrValue("ID") + " not found");
    }
    return parent;
}

/*
 * Process a FF3 record for a gene, which uses the explicit tree.
 * Return the new leaf node.
 */
GxfFeatureNode* GxfFeatureTree::loadGff3GeneRecord(const GxfFeature* gxfFeature,
                                                   GxfFeatureNode* geneTreeLeaf) const {
    GxfFeatureNode* parent = findGff3Parent(geneTreeLeaf, gxfFeature);
    GxfFeatureNode* child = new GxfFeatureNode(gxfFeature);
    parent->addChild(child);
    return child;
}

/* Get the desired type the parent of GTF feature
 * WARNING: this assumed that the hierarchy is
 * gene->transcript->{everything else}
 * FIXME: this is what GFF3 does, which might not be right.
 */
const string& GxfFeatureTree::getGtfParentType(const string& featureType) const {
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
GxfFeatureNode* GxfFeatureTree::findGtfParent(GxfFeatureNode* geneTreeLeaf,
                                              const GxfFeature* gxfFeature) const {
    const string& parentType = getGtfParentType(gxfFeature->fType);
    GxfFeatureNode* parent = geneTreeLeaf;
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
RemapStatus GxfFeatureNode::calcRemapStatus(bool srcSeqInMapping) const {
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
        if (fMappedFeatures.size() == 1) {
            return REMAP_STATUS_PARTIAL_CONTIG;
        } else {
            return REMAP_STATUS_PARTIAL_FRAGMENT;
        }
    }
}

/* print node for debugging */
void GxfFeatureNode::dumpNode(ostream& fh) const {
    const string status = remapStatusToStr(fRemapStatus);
    fh << "src" << "\t" << ((fFeature == NULL) ? "NULL" : fFeature->toString()) << endl;
    for (int i = 0; i < fAllFeatures.size(); i++) {
        fh << (fMappedFeatures.contains(fAllFeatures[i]) ? "mapped" : "unmapped") << "\t" << status << "\t" << fAllFeatures[i]->toString() << endl;
    }
}

/* recursively print for debugging */
void GxfFeatureNode::dump(ostream& fh) const {
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
GxfFeatureNode* GxfFeatureTree::loadGtfGeneRecord(const GxfFeature* gxfFeature,
                                                  GxfFeatureNode* geneTreeLeaf) const {
    GxfFeatureNode* parent = findGtfParent(geneTreeLeaf, gxfFeature);
    GxfFeatureNode* child = new GxfFeatureNode(gxfFeature);
    parent->addChild(child);
    return child;
}

/*
 * Process a GxfRecord for a gene, return False if no more for this gene.
 */
bool GxfFeatureTree::loadGeneRecord(GxfParser *gxfParser,
                                    const GxfRecord* gxfRecord,
                                    GxfFeatureNode* geneTreeRoot,
                                    GxfFeatureNode*& geneTreeLeaf,
                                    GxfRecordVector& queuedRecords) const {
    if (instanceOf(gxfRecord, GxfLine)) {
        queuedRecords.push_back(gxfRecord);
        return true;
    } else {
        const GxfFeature* gxfFeature = dynamic_cast<const GxfFeature*>(gxfRecord);
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
 * of genes to be moved to the end, but GENCODE doesn't do this.  This whole
 * thing is annoying due to the lack of explicit structure in GTF.
 */
GxfFeatureNode* GxfFeatureTree::loadGene(GxfParser *gxfParser,
                                         const GxfFeature* geneFeature) {
    assert(geneFeature->fType == GxfFeature::GENE);

    GxfFeatureNode* geneTreeRoot = new GxfFeatureNode(geneFeature);
    GxfFeatureNode* geneTreeLeaf = geneTreeRoot;  // were we are currently working
    GxfRecordVector queuedRecords;
    const GxfRecord* gxfRecord = NULL;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (not loadGeneRecord(gxfParser, gxfRecord, geneTreeRoot, geneTreeLeaf, queuedRecords)) {
            break;
        }
    }
    queueRecords(gxfParser, queuedRecords);
    return geneTreeRoot;
}

/* constructor */
GxfFeatureTree::GxfFeatureTree(GxfParser *gxfParser,
                               const GxfFeature* geneFeature):
    fGene(NULL) {
    fGene = loadGene(gxfParser, geneFeature);
}

/* Destructor */
GxfFeatureTree::~GxfFeatureTree() {
    delete fGene;
}

/* print for debugging */
void GxfFeatureTree::dump(ostream& fh) const {
    fGene->dump(fh);
}
