/*
 * Tree structure use to store genes
 */
#include "featureTree.hh"
#include <iostream>
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

/* ensembl non-coding gene biotypes to skip */
static const char* automaticNonCodingGeneBiotypes[] = {
    "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ribozyme", "rRNA", "scaRNA",
    "snoRNA", "snRNA", "sRNA", NULL
};

/* is ensembl small non-coding gene */
bool FeatureNode::isAutomaticSmallNonCodingGene() const {
    if (fFeature->fSource != "ENSEMBL") {
        return false;
    }
    const string& bioType = fFeature->getTypeBiotype();
    for (int i = 0; automaticNonCodingGeneBiotypes[i] != NULL; i++) {
        if (bioType == automaticNonCodingGeneBiotypes[i]) {
            return true;
        }
    }
    return false;
}

/* set the remap number of mappings attribute on this node.  not recursive,
 * since it's only set on gene/transcript */
void FeatureNode::setNumMappingsAttr() {
    fFeature->getAttrs().update(AttrVal(REMAP_NUM_MAPPINGS_ATTR, toString(fNumMappings)));
}

/* recursively set the remap status attribute */
void FeatureNode::rsetRemapStatusAttr() {
    if (fRemapStatus != REMAP_STATUS_NONE) {
        fFeature->getAttrs().update(AttrVal(REMAP_STATUS_ATTR, remapStatusToStr(fRemapStatus)));
    }
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetRemapStatusAttr();
    }
}

/* recursively set the target status attribute */
void FeatureNode::rsetTargetStatusAttr() {
    if (isGeneOrTranscript()) {
        if (fTargetStatus != TARGET_STATUS_NA) {
            fFeature->getAttrs().update(AttrVal(REMAP_TARGET_STATUS_ATTR, targetStatusToStr(fTargetStatus)));
        }
        for (int i = 0; i < fChildren.size(); i++) {
            fChildren[i]->rsetTargetStatusAttr();
        }
    }
}

/* recursively set the target status attribute on fFeature node. */
void FeatureNode::rsetSubstitutedMissingTargetAttr(const string& targetVersion) {
    fFeature->getAttrs().update(AttrVal(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR, targetVersion));
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetSubstitutedMissingTargetAttr(targetVersion);
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
                                      const GxfFeature* feature) {
    const string& parentId = feature->getAttrValue(GxfFeature::PARENT_ATTR);
    FeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->getAttrValue(GxfFeature::ID_ATTR) != parentId)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node " +  parentId + " for " + feature->getAttrValue(GxfFeature::ID_ATTR) + " not found");
    }
    return parent;
}

/*
 * Process a FF3 record for a gene, which uses the explicit tree.
 * Return the new leaf node.
 */
FeatureNode* GeneTree::loadGff3GeneRecord(GxfFeature* feature,
                                          FeatureNode* geneTreeLeaf) {
    FeatureNode* parent = findGff3Parent(geneTreeLeaf, feature);
    FeatureNode* child = new FeatureNode(feature);
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
                                     const GxfFeature* feature) {
    const string& parentType = getGtfParentType(feature->fType);
    FeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->fType != parentType)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node of type " + parentType + "  not found for type " + parent->fFeature->fType);
    }
    return parent;
}

/* compute status for a transmapped feature.  this only looks at a single
 * level of mappings, not a tree. srcSeqInMapping indicates of the srcSequence
 * was in the genomic map */
RemapStatus TransMappedFeature::calcRemapStatus(bool srcSeqInMapping) const {
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

/* recursively set the remap status. */
void FeatureNode::rsetRemapStatus(RemapStatus remapStatus) {
    fRemapStatus = remapStatus;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetRemapStatus(remapStatus);
    }
}

/* Set the target status. Not recursive */
void FeatureNode::setTargetStatus(TargetStatus targetStatus) {
    fTargetStatus = targetStatus;
}

/* Recursively set the target status. */
void FeatureNode::rsetTargetStatus(TargetStatus targetStatus) {
    fTargetStatus = targetStatus;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetTargetStatus(targetStatus);
    }
}

/* do any child belond to the specified status */
bool FeatureNode::anyChildWithRemapStatus(unsigned remapStatusSet) const {
    for (size_t i = 0; i < fChildren.size(); i++) {
        if ((fChildren[i]->fRemapStatus & remapStatusSet) != 0) {
            return true;
        }
    }
    return false;
}

/* do all child have belong to the specified status set */
bool FeatureNode::allChildWithRemapStatus(unsigned remapStatusSet) const {
    for (size_t i = 0; i < fChildren.size(); i++) {
        if ((fChildren[i]->fRemapStatus & remapStatusSet) == 0) {
            return false;
        }
    }
    return true;
}

/* get the size of a transcript, in exons */
int FeatureNode::getTranscriptExonSize() const {
    assert(isTranscript());
    int size = 0;
    for (int iExon = 0; iExon < fChildren.size(); iExon++) {
        if (fChildren[iExon]->isExon()) {
            size += fChildren[iExon]->fFeature->size();
        }
    }
    return size;
}

/* count overlapping bases */
int FeatureNode::getOverlapAmount(const FeatureNode* other) const {
    int maxStart = max(other->fFeature->fStart, fFeature->fStart);
    int minEnd = min(other->fFeature->fEnd, fFeature->fEnd);
    return (maxStart <= minEnd) ? (minEnd-maxStart)+1 : 0;
}

/* count exon overlap with exons of another transcript */
int FeatureNode::countExonOverlap(const FeatureNode* exon1,
                                  const FeatureNode* trans2) const {
    int totalOverlap = 0;
    for (int iFeat2 = 0; iFeat2 < trans2->fChildren.size(); iFeat2++) {
        if (trans2->fChildren[iFeat2]->isExon()) {
            const FeatureNode* exon2 = trans2->fChildren[iFeat2];
            totalOverlap += exon1->getOverlapAmount(exon2);
        }
    }
    return totalOverlap;
}

/* get exon similarity */
float FeatureNode::getExonSimilarity(const FeatureNode* trans2) const {
    assert(isTranscript());
    assert(trans2->isTranscript());
    int totalOverlap = 0;
    for (int iFeat1 = 0; iFeat1 < fChildren.size(); iFeat1++) {
        if (fChildren[iFeat1]->isExon()) {
            totalOverlap += countExonOverlap(fChildren[iFeat1], trans2);
        }
    }
    return float(2*totalOverlap)/float(getTranscriptExonSize() + trans2->getTranscriptExonSize());
}

/* get the maximum transcript similarity for a gene */
float FeatureNode::getMaxTranscriptSimilarity(const FeatureNode* gene2) const {
    assert(isGene());
    assert(gene2->isGene());
    float maxSimilarity = 0.0;
    for (int iTrans1 = 0; (iTrans1 < fChildren.size()) && (maxSimilarity < 1.0); iTrans1++) {
        for (int iTrans2 = 0; (iTrans2 < gene2->fChildren.size()) && (maxSimilarity < 1.0); iTrans2++) {
            float similarity = fChildren[iTrans1]->getExonSimilarity(gene2->fChildren[iTrans2]);
            maxSimilarity = max(similarity, maxSimilarity);
        }
    }
    return maxSimilarity;
}

/* do any child belond to the specified status */
bool ResultFeatureTrees::anyChildWithRemapStatus(unsigned remapStatusSet) const {
    return ((mapped != NULL) && mapped->anyChildWithRemapStatus(remapStatusSet))
        || ((unmapped != NULL) && unmapped->anyChildWithRemapStatus(remapStatusSet));
}

/* do all child have belong to the specified status set */
bool ResultFeatureTrees::allChildWithRemapStatus(unsigned remapStatusSet) const {
    return ((mapped == NULL) || ((mapped != NULL) && mapped->allChildWithRemapStatus(remapStatusSet)))
        && ((unmapped == NULL) || ((unmapped != NULL) && unmapped->allChildWithRemapStatus(remapStatusSet)));
}

/* determine gene or transcript remap status from children.  This doesn't
 * handle GENE_CONFLICT or_GENE_SIZE_CHANGE, which are forced.
 */
RemapStatus ResultFeatureTrees::calcBoundingFeatureRemapStatus(bool srcSeqInMapping) const {
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

/* clone tree */
FeatureNode* FeatureNode::clone() const {
    FeatureNode *newNode = new FeatureNode(fFeature->clone());
    for (int i = 0; i < fChildren.size(); i++) {
        newNode->fChildren.push_back(fChildren[i]->clone());
    }
    return newNode;
}

/* print node for debugging */
void FeatureNode::dumpNode(ostream& fh) const {
    fh << fFeature->toString() << endl;
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
        GxfFeature* feature = dynamic_cast<GxfFeature*>(gxfRecord);
        if (feature->fType == GxfFeature::GENE) {
            queuedRecords.push_back(gxfRecord); // next gene
            return false;
        } else {
            // FIXME: should have parser/format convert to common in memory object
            if (gxfParser->getFormat() == GFF3_FORMAT) {
                geneTreeLeaf = loadGff3GeneRecord(feature, geneTreeLeaf);
            } else {
                geneTreeLeaf = loadGtfGeneRecord(feature, geneTreeLeaf);
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

