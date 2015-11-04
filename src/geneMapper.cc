#include "geneMapper.hh"
#include "featureMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "transcriptMapper.hh"
#include "targetAnnotations.hh"

/* fraction of gene expansion that causes a rejection */
const float geneExpansionThreshold = 0.50;

/*  mapinfo TSV headers, terminated by NULL */
static const char* mappingInfoHeaders[] = {
    "id", "name", "type", "biotype", "source",
    "srcChrom", "srcStart", "srcEnd", "srcStrand",
    "mappedChrom", "mappedStart", "mappedEnd", "mappedStrand",
    "mappingStatus", "numMappings",
    "targetStatus", "targetBiotype", "targetSubst", NULL
};

/* ensembl non-coding gene biotypes to skip */
static const char* automaticNonCodingGeneBiotypes[] = {
    "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ribozyme", "rRNA", "scaRNA",
    "snoRNA", "snRNA", "sRNA", NULL
};

/* is ensembl small non-coding gene */
bool GeneMapper::isAutomaticSmallNonCodingGene(FeatureNode* geneTree) {
    if (geneTree->fFeature->fSource != "ENSEMBL") {
        return false;
    }
    const string& bioType = geneTree->fFeature->getTypeBiotype();
    for (int i = 0; automaticNonCodingGeneBiotypes[i] != NULL; i++) {
        if (bioType == automaticNonCodingGeneBiotypes[i]) {
            return true;
        }
    }
    return false;
}
/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const GxfFeature* feature) const {
    return fGenomeTransMap->haveQuerySeq(feature->fSeqid);
}

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(FeatureNode* featureNode) const {
    return isSrcSeqInMapping(featureNode->fFeature);
}

/* process one transcript */
void GeneMapper::processTranscript(FeatureNode* transcriptTree,
                                   ostream* transcriptPslFh) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcriptTree, fTargetAnnotations,
                                      isSrcSeqInMapping(transcriptTree), transcriptPslFh);
    transcriptMapper.mapTranscriptFeatures(transcriptTree);
    transcriptTree->fTargetStatus = getTargetAnnotationStatus(transcriptTree);
    if (transcriptTree->fTargetStatus == TARGET_STATUS_NONOVERLAP) {
        forceToUnmappedDueToTargetStatus(transcriptTree, transcriptTree->fTargetStatus);
    }
}

/* process all transcripts of gene. */
void GeneMapper::processTranscripts(FeatureNode* geneTree,
                                    ostream* transcriptPslFh) const {
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        processTranscript(transcriptTree, transcriptPslFh);
    }
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when we discover conflicts at the gene level. */
void GeneMapper::forceToUnmappedDueToRemapStatus(FeatureNode* featureNode,
                                                 RemapStatus remapStatus) const {
    FeatureMapper::forceToUnmapped(featureNode);
    featureNode->setRemapStatus(remapStatus);
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        forceToUnmappedDueToRemapStatus(featureNode->fChildren[i], remapStatus);
    }
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when target status indicates a bad mapping */
void GeneMapper::forceToUnmappedDueToTargetStatus(FeatureNode* featureNode,
                                                  TargetStatus targetStatus) const {
    FeatureMapper::forceToUnmapped(featureNode);
    featureNode->fTargetStatus = targetStatus;
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        forceToUnmappedDueToTargetStatus(featureNode->fChildren[i], targetStatus);
    }
}

/* are there any mapped transcripts? */
bool GeneMapper::haveMappedTranscripts(FeatureNode* geneTree) const {
    return (geneTree->fChildren.size() > 0) and (geneTree->fChildren[0]->fMappedFeatures.size() > 0);
}

/* are there any unmapped transcripts? */
bool GeneMapper::haveUnmappedTranscripts(FeatureNode* geneTree) const {
    return (geneTree->fChildren.size() > 0) and (geneTree->fChildren[0]->fUnmappedFeatures.size() > 0);
}

/* check if gene contains transcripts with different mapped seqid/strand  */
bool GeneMapper::hasMixedMappedSeqStrand(FeatureNode* geneTree) const {
    string seqid, strand;
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fMappedFeatures.size() > 0) {
            GxfFeature* mappedTranscript = transcriptTree->fMappedFeatures[0];
            if (seqid == "") {
                seqid = mappedTranscript->fSeqid;
                strand = mappedTranscript->fStrand;
            } else if ((mappedTranscript->fSeqid != seqid)
                       or (mappedTranscript->fStrand != strand)) {
                return true;
            }
        }
    }
    return false;
}

/* check if gene contains transcripts of a gene have a target status of nonoverlap */
bool GeneMapper::hasTargetStatusNonOverlap(FeatureNode* geneTree) const {
    if (geneTree->fTargetStatus == TARGET_STATUS_NONOVERLAP) {
        // handles weird case of ENSG00000239810 were gene doesn't overlap
        // however the only transcript is `new'.
        return true;
    }
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        if (geneTree->fChildren[i]->fTargetStatus == TARGET_STATUS_NONOVERLAP) {
            return true;
        }
    }
    return false;
}

/* compute the length of a mapped gene */
int GeneMapper::calcMappedGeneLength(FeatureNode* geneTree) const {
    int start = 0, end = 0;
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fMappedFeatures.size() > 0) {
            GxfFeature* mappedTranscript = transcriptTree->fMappedFeatures[0];
            if (start == 0) {
                start = mappedTranscript->fStart;
                end = mappedTranscript->fEnd;
            } else {
                start = min(start, mappedTranscript->fStart);
                end = max(end, mappedTranscript->fEnd);
            }
        }
    }
    return (end == 0) ? 0 : (end-start)+1;
}

/* check if gene transcripts cause the gene size to expand beyond a
 * threshold*/
bool GeneMapper::hasExcessiveSizeChange(FeatureNode* geneTree) const {
    int srcGeneLength = geneTree->fFeature->size();
    int mappedGeneLength = calcMappedGeneLength(geneTree);
    if (mappedGeneLength == 0) {
        return 0;  // not mapped
    } else {
        return (abs(float(mappedGeneLength-srcGeneLength)) / float(srcGeneLength)) > geneExpansionThreshold;
    }
}

/* set the number of mappings of the gene  */
void GeneMapper::setNumGeneMappings(FeatureNode* geneTree) const {
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        geneTree->fNumMappings = max(geneTree->fNumMappings, geneTree->fChildren[i]->fNumMappings);
    }
}

/* update gene bounds given a mapped transcript */
void GeneMapper::updateMappedGeneBounds(FeatureNode* transcriptTree,
                                        string& seqid, string& strand,
                                        int& start, int& end) const {
    assert(transcriptTree->fMappedFeatures.size() == 1);
    GxfFeature* mappedTranscript = transcriptTree->fMappedFeatures[0];
    assert(mappedTranscript->fType == GxfFeature::TRANSCRIPT);
    if (seqid == "") {
        // first
        seqid = mappedTranscript->fSeqid;
        strand = mappedTranscript->fStrand;
        start = mappedTranscript->fStart;
        end = mappedTranscript->fEnd;
    } else {
        start = min(start, mappedTranscript->fStart);
        end = max(end, mappedTranscript->fEnd);
    }
}


/* should we substitute target version of gene?  */
bool GeneMapper::shouldSubstituteMissingTarget(FeatureNode* geneTree) const {
    // requested and a mis-mapped gene or the right biotype to a mapped sequence
    if (not ((fSubstituteMissingTargetVersion.size() > 0)
             and ((geneTree->fTargetStatus == TARGET_STATUS_NONOVERLAP)
                  or (geneTree->fTargetStatus == TARGET_STATUS_LOST)))) {
        return false; // not substituting or not right target status
    }
    const FeatureNode* targetGene = getTargetAnnotationNode(geneTree);
    if (targetGene == NULL) {
        // this should not have happened, but it does because of a transcript id
        // ENST00000426406 incorrectly being moved to a different gene
        return false;
    }
    if (not isSrcSeqInMapping(targetGene->fFeature)) {
        return false;  // sequence not being mapped (moved chroms)
    }

    return (targetGene->fFeature->getTypeBiotype() == geneTree->fFeature->getTypeBiotype());
}

/* copy target gene to substitute for this mapping and set attributes  */
void GeneMapper::substituteMissingTarget(FeatureNode* geneTree) const {
    assert(not haveMappedTranscripts(geneTree));
    const FeatureNode* targetGene = getTargetAnnotationNode(geneTree);
    geneTree->fSubstitutedMissingTarget = targetGene->clone(geneTree->fFeature->getFormat());
    // Set flag in both trees
    geneTree->fSubstitutedMissingTarget->setSubstitutedMissingTargetAttrOnFeature(fSubstituteMissingTargetVersion);
    geneTree->setSubstitutedMissingTargetAttrOnUnmapped(fSubstituteMissingTargetVersion);
}

/* If there are any mapped transcripts of a gene, add a gene
 * record for the bounds */
void GeneMapper::buildMappedGeneFeature(FeatureNode* geneTree,
                                        bool srcSeqInMapping) const {
    // calculate bounds
    string seqid, strand;
    int start = 0, end = 0;
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fMappedFeatures.size() > 0) {
            updateMappedGeneBounds(transcriptTree, seqid, strand, start, end);
        }
    }
    FeatureMapper::mapBounding(geneTree, seqid, start-1, end, strand);  // takes zero bases
}

/* If there are any unmapped transcripts of a gene, add a gene
 * record for the bounds */
void GeneMapper::buildUnmappedGeneFeature(FeatureNode* geneTree,
                                          bool srcSeqInMapping) const {
    FeatureMapper::mapBounding(geneTree);
}

/* Build gene features */
void GeneMapper::buildGeneFeature(FeatureNode* geneTree) const {
    bool srcSeqInMapping = isSrcSeqInMapping(geneTree);
    setNumGeneMappings(geneTree);
    if (haveMappedTranscripts(geneTree)) {
        buildMappedGeneFeature(geneTree, srcSeqInMapping);
    }
    if (haveUnmappedTranscripts(geneTree)) {
        buildUnmappedGeneFeature(geneTree, srcSeqInMapping);
    }
    geneTree->fTargetStatus = getTargetAnnotationStatus(geneTree);
    if ((geneTree->fTargetStatus == TARGET_STATUS_LOST)
        && hasTargetStatusNonOverlap(geneTree)) {
        // lost because at least one transcript doesn't overlap the target
        geneTree->fTargetStatus = TARGET_STATUS_NONOVERLAP;
    }
}

/* output GFF3 mapped ##sequence-region if not already written */
void GeneMapper::outputMappedSeqRegionIfNeed(FeatureNode* geneTree,
                                             ostream& mappedGxfFh) {
    assert((geneTree->fSubstitutedMissingTarget != NULL) || (geneTree->fMappedFeatures.size() > 0));
    if (geneTree->fFeature->getFormat() == GFF3_FORMAT) {
        // count be due to substituted gene
        const string& mappedSeqId = (geneTree->fSubstitutedMissingTarget != NULL)
            ? geneTree->fSubstitutedMissingTarget->fFeature->fSeqid
            : geneTree->fMappedFeatures[0]->fSeqid;
        if (not checkRecordSeqRegionWritten(mappedSeqId)) {
            mappedGxfFh << "##sequence-region " << mappedSeqId << " 1 "
                        << fGenomeTransMap->getTargetSeqSize(mappedSeqId) << endl;
        }
    }
}

/*
 * recursive output of a GxF mapped feature tree
 */
void GeneMapper::outputMapped(FeatureNode* featureNode,
                              ostream& mappedGxfFh) const {
    for (int i = 0; i < featureNode->fMappedFeatures.size(); i++) {
        mappedGxfFh << featureNode->fMappedFeatures[i]->toString() << endl;
    }
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        outputMapped(featureNode->fChildren[i], mappedGxfFh);
    }
}

/*
 * recursive output of a GxF unmapped feature tree
 */
void GeneMapper::outputUnmapped(FeatureNode* featureNode,
                                ostream& unmappedGxfFh) const {
    for (int i = 0; i < featureNode->fUnmappedFeatures.size(); i++) {
        unmappedGxfFh << featureNode->fUnmappedFeatures[i]->toString() << endl;
    }
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        outputUnmapped(featureNode->fChildren[i], unmappedGxfFh);
    }
}

/*
 * recursive output of a GxF substituted feature tree
 */
void GeneMapper::outputSubstituted(FeatureNode* featureNode,
                                   ostream& mappedGxfFh) const {
    mappedGxfFh << featureNode->fFeature->toString() << endl;
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        outputSubstituted(featureNode->fChildren[i], mappedGxfFh);
    }
}

/* output genes */
void GeneMapper::outputFeatures(FeatureNode* geneTree,
                                ostream& mappedGxfFh,
                                ostream& unmappedGxfFh) {
    if (haveMappedTranscripts(geneTree)) {
        outputMappedSeqRegionIfNeed(geneTree, mappedGxfFh);
        outputMapped(geneTree, mappedGxfFh);
    }
    if (haveUnmappedTranscripts(geneTree)) {
        outputUnmapped(geneTree, unmappedGxfFh);
    }
    if (geneTree->fSubstitutedMissingTarget != NULL) {
        outputMappedSeqRegionIfNeed(geneTree, mappedGxfFh);
        outputSubstituted(geneTree->fSubstitutedMissingTarget, mappedGxfFh);
    }
}

/* output info TSV header */
void GeneMapper::outputInfoHeader(ostream& mappingInfoFh) const {
    for (int i = 0; mappingInfoHeaders[i] != NULL; i++) {
        if (i > 0) {
            mappingInfoFh << "\t";
        }
        mappingInfoFh << mappingInfoHeaders[i];
    }
    mappingInfoFh << endl;
}

/* get target annotation for a feature, if available */
const GxfFeature* GeneMapper::getTargetAnnotation(FeatureNode* featureNode) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    } else {
        return fTargetAnnotations->getFeature(featureNode->fFeature->getTypeId(),
                                              featureNode->fFeature->fSeqid);
    }
}

/* get target annotation node for a feature, if available */
const FeatureNode* GeneMapper::getTargetAnnotationNode(FeatureNode* featureNode) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    } else {
        return fTargetAnnotations->getFeatureNode(featureNode->fFeature->getTypeId(),
                                                  featureNode->fFeature->fSeqid);
    }
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
TargetStatus GeneMapper::getTargetAnnotationStatus(FeatureNode* featureNode) const {
    const GxfFeature* targetFeature = getTargetAnnotation(featureNode);
    if (fTargetAnnotations == NULL) {
        return TARGET_STATUS_NA;
    }
    if (targetFeature == NULL) {
        return TARGET_STATUS_NEW;
    }
    if (featureNode->fMappedFeatures.size() == 0) {
        return TARGET_STATUS_LOST;
    }
    assert(featureNode->fMappedFeatures.size() == 1);  // genes/transcripts don't split
    if (featureNode->fMappedFeatures[0]->overlaps(targetFeature)) {
        return TARGET_STATUS_OVERLAP;
    } else {
        return TARGET_STATUS_NONOVERLAP;
    }
}

/* If target gene annotations are available, get biotype of target feature */
const string& GeneMapper::getTargetAnnotationBiotype(FeatureNode* featureNode) const {
    const GxfFeature* targetFeature = getTargetAnnotation(featureNode);
    if (targetFeature == NULL) {
        return emptyString;
    } else {
        return targetFeature->getTypeBiotype();
    }
}

/* output info on one bounding feature */
void GeneMapper::outputFeatureInfo(FeatureNode* featureNode,
                                   bool substituteMissingTarget,
                                   ostream& mappingInfoFh) const {
    assert(featureNode->fMappedFeatures.size() <= 1);  // only handles bounding features
    mappingInfoFh << featureNode->fFeature->getTypeId() << "\t"
                  << featureNode->fFeature->getTypeName() << "\t"
                  << featureNode->fFeature->fType << "\t"
                  << featureNode->fFeature->getTypeBiotype() << "\t"
                  << featureNode->fFeature->fSource << "\t"
                  << featureNode->fFeature->fSeqid << "\t"
                  << featureNode->fFeature->fStart << "\t"
                  << featureNode->fFeature->fEnd << "\t"
                  << featureNode->fFeature->fStrand << "\t";
    if (featureNode->fMappedFeatures.size() > 0) {
        const GxfFeature* mappedFeature = featureNode->fMappedFeatures[0];
        mappingInfoFh << mappedFeature->fSeqid << "\t"
                      << mappedFeature->fStart << "\t"
                      << mappedFeature->fEnd << "\t"
                      << mappedFeature->fStrand << "\t";
    } else {
        mappingInfoFh << "\t0\t0\t.\t";
    }
    mappingInfoFh << remapStatusToStr(featureNode->fRemapStatus) << "\t"
                  << featureNode->fNumMappings << "\t"
                  << targetStatusToStr(featureNode->fTargetStatus) << "\t"
                  << getTargetAnnotationBiotype(featureNode) << "\t"
                  << (substituteMissingTarget ? fSubstituteMissingTargetVersion : emptyString)
                  << endl;
}

/*
 * Output information about gene mapping
 */
void GeneMapper::outputInfo(FeatureNode* geneNode,
                            ostream& mappingInfoFh) const {
    bool substituteMissingTarget = geneNode->fSubstitutedMissingTarget != NULL;
    outputFeatureInfo(geneNode, substituteMissingTarget, mappingInfoFh);
    // transcripts
    for (int i = 0; i < geneNode->fChildren.size(); i++) {
        outputFeatureInfo(geneNode->fChildren[i], substituteMissingTarget, mappingInfoFh);
    }
}

/*
 * map and output one gene's annotations
 */
void GeneMapper::processGene(GxfParser *gxfParser,
                             GxfFeature* geneFeature,
                             ostream& mappedGxfFh,
                             ostream& unmappedGxfFh,
                             ostream& mappingInfoFh,
                             ostream* transcriptPslFh) {
    FeatureNode* geneTree = GeneTree::geneTreeFactory(gxfParser, geneFeature);
    processTranscripts(geneTree, transcriptPslFh);
    // handle gene level issues
    if (fSkipAutomaticNonCoding and 
        isAutomaticSmallNonCodingGene(geneTree)) {
        forceToUnmappedDueToRemapStatus(geneTree, REMAP_STATUS_AUTOMATIC_NON_CODING);
    }
    if (hasMixedMappedSeqStrand(geneTree)) {
        forceToUnmappedDueToRemapStatus(geneTree, REMAP_STATUS_GENE_CONFLICT);
    } else if (hasExcessiveSizeChange(geneTree)) {
        forceToUnmappedDueToRemapStatus(geneTree, REMAP_STATUS_GENE_SIZE_CHANGE);
    }
    buildGeneFeature(geneTree);
    if (hasTargetStatusNonOverlap(geneTree)) {
        // must come after building gene do to weird case of ENSG00000239810
        // were gene doesn't overlap however the only transcript is `new'.
        forceToUnmappedDueToTargetStatus(geneTree, TARGET_STATUS_NONOVERLAP);
    }
    geneTree->setBoundingFeatureRemapStatus();
    geneTree->setRemapStatusAttr();
    geneTree->setNumMappingsAttr();
    geneTree->setTargetStatusAttr();

    // must be done after forcing status above
    if (shouldSubstituteMissingTarget(geneTree)) {
        substituteMissingTarget(geneTree);
    }
    outputFeatures(geneTree, mappedGxfFh, unmappedGxfFh);
    outputInfo(geneTree, mappingInfoFh);
    delete geneTree;
}

/* process a record, this may consume additional feature records  */
void GeneMapper::processRecord(GxfParser *gxfParser,
                               GxfRecord* gxfRecord,
                               ostream& mappedGxfFh,
                               ostream& unmappedGxfFh,
                               ostream& mappingInfoFh,
                               ostream* transcriptPslFh) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        processGene(gxfParser, dynamic_cast<GxfFeature*>(gxfRecord), mappedGxfFh, unmappedGxfFh, mappingInfoFh, transcriptPslFh);
    } else if ((gxfParser->getGxfFormat() == GFF3_FORMAT) and (gxfRecord->toString().find("##sequence-region") == 0)) {
        // mapped ##sequence-region records are created based on sequences written
        unmappedGxfFh << gxfRecord->toString() << endl;
        delete gxfRecord;
    } else {
        // write all others as-is
        mappedGxfFh << gxfRecord->toString() << endl;
        unmappedGxfFh << gxfRecord->toString() << endl;
        delete gxfRecord;
    }
}

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfParser *gxfParser,
                        ostream& mappedGxfFh,
                        ostream& unmappedGxfFh,
                        ostream& mappingInfoFh,
                        ostream* transcriptPslFh) {
    outputInfoHeader(mappingInfoFh);
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        processRecord(gxfParser, gxfRecord, mappedGxfFh, unmappedGxfFh, mappingInfoFh, transcriptPslFh);
    }
}

