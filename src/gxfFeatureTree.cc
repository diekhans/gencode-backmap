/*
 * Tree structure use to store genes
 */
#include "gxfFeatureTree.hh"
#include <ostream>

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
    GxfFeatureNode* child = new GxfFeatureNode(gxfFeature, REMAP_STATUS_NONE);
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

/*
 * Process a GTF record for a gene, which uses knowledge of
 * the GENCODE structure to reproduce the hierarchy.
 * Return the new leaf node.
 */
GxfFeatureNode* GxfFeatureTree::loadGtfGeneRecord(const GxfFeature* gxfFeature,
                                                  GxfFeatureNode* geneTreeLeaf) const {
    GxfFeatureNode* parent = findGtfParent(geneTreeLeaf, gxfFeature);
    GxfFeatureNode* child = new GxfFeatureNode(gxfFeature, REMAP_STATUS_NONE);
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

    GxfFeatureNode* geneTreeRoot = new GxfFeatureNode(geneFeature, REMAP_STATUS_NONE);
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

/* Attribute name used for remap status */
const string REMAP_ATTR_NAME = "remap_status";

/* convert a remap status to a AttrVal  */
const AttrVal& remapStatusToAttrVal(RemapStatus remapStatus) {
    switch (remapStatus) {
    case REMAP_STATUS_NONE: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_none");
        return status;
    }
    case REMAP_STATUS_FULL_CONTIG: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_full_contig");
        return status;
    }
    case REMAP_STATUS_FULL_FRAGMENT: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_full_fragment");
        return status;
    }
    case REMAP_STATUS_PARTIAL_CONTIG: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_partial_contig");
        return status;
    }
    case REMAP_STATUS_PARTIAL_FRAGMENT: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_partial_fragment");
        return status;
    }
    case REMAP_STATUS_DELETED: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_deleted");
        return status;
    }
    case REMAP_STATUS_NO_SEQ_MAP: {
        static const AttrVal status(REMAP_ATTR_NAME, "remap_status_no_seq_map");
        return status;
    }
    default:
        throw invalid_argument("BUG: invalid remap status");
    }
}
