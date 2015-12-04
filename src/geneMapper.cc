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
bool GeneMapper::isAutomaticSmallNonCodingGene(const FeatureNode* geneTree) {
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
bool GeneMapper::isSrcSeqInMapping(const FeatureNode* featureNode) const {
    return isSrcSeqInMapping(featureNode->fFeature);
}

/* process one transcript */
ResultFeatureTrees GeneMapper::processTranscript(const FeatureNode* transcriptTree,
                                                 ostream* transcriptPslFh) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcriptTree, fTargetAnnotations,
                                      isSrcSeqInMapping(transcriptTree), transcriptPslFh);
    ResultFeatureTrees mappedTranscript = transcriptMapper.mapTranscriptFeatures(transcriptTree);
    TargetStatus targetStatus = getTargetAnnotationStatus(&mappedTranscript);
    mappedTranscript.setTargetStatus(targetStatus);
    return mappedTranscript;
}

/* process all transcripts of gene. */
ResultFeatureTreesVector GeneMapper::processTranscripts(const FeatureNode* geneTree,
                                                        ostream* transcriptPslFh) const {
    ResultFeatureTreesVector mappedTranscripts;
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        mappedTranscripts.push_back(processTranscript(transcriptTree, transcriptPslFh));
    }
    return mappedTranscripts;
}

/* find a matching gene or transcript node given node by id */
FeatureNode* GeneMapper::findMatchingBoundingNode(const FeatureNodeVector& features,
                                                  const FeatureNode* feature) const {
    assert(feature->isGeneOrTranscript());
    for (int i = 0; i < features.size(); i++)  {
        assert(features[i]->isGeneOrTranscript());
        if (features[i]->fFeature->getTypeId() == feature->fFeature->getTypeId()) {
            return features[i];  // found match!
        }
    }
    return NULL;
}

/* copy mapping metadata */
void GeneMapper::copyMappingMetadata(const FeatureNode* origFeature,
                                     FeatureNode* newFeature) const {
    newFeature->setRemapStatus(origFeature->fRemapStatus);
    newFeature->setTargetStatus(origFeature->fTargetStatus);
    newFeature->fNumMappings = max(origFeature->fNumMappings, newFeature->fNumMappings);
}

/* copy gene and transcript attributes to a fresh unmapped clone,
 * only at the transcript and gene levels. */
void GeneMapper::copyGeneMetadata(const FeatureNode* origGene,
                                  FeatureNode* newGene) const {
    copyMappingMetadata(origGene, newGene);

    // copy transcript metadata
    for (size_t i = 0; i < newGene->fChildren.size(); i++) {
        FeatureNode* newTranscript = newGene->fChildren[i];
        const FeatureNode* origTranscript = findMatchingBoundingNode(origGene->fChildren, newTranscript);
        if (origTranscript != NULL) {
            copyMappingMetadata(origTranscript, newTranscript);
        }
    }
}

/* force to unmapped */
void GeneMapper::forceToUnmapped(ResultFeatureTrees* mappedGene) const {
    // create copy and update attributes at gene/transcript level before freeing.
    FeatureNode* srcCopy = mappedGene->src->clone(mappedGene->src->fFeature->getFormat());
    if (mappedGene->mapped != NULL) {
        copyGeneMetadata(mappedGene->mapped, srcCopy);
    }
    if (mappedGene->unmapped != NULL) {
        copyGeneMetadata(mappedGene->unmapped, srcCopy);
    }
    srcCopy->rsetRemapStatusAttr();
    srcCopy->setNumMappingsAttr();
    srcCopy->rsetTargetStatusAttr();
    mappedGene->freeMapped();
    mappedGene->freeUnmapped();
    mappedGene->unmapped = srcCopy;
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when we discover conflicts at the gene level. */
void GeneMapper::forceToUnmappedDueToRemapStatus(ResultFeatureTrees* mappedGene,
                                                 RemapStatus remapStatus) const {
    forceToUnmapped(mappedGene);
    mappedGene->rsetRemapStatus(remapStatus);
    mappedGene->rsetRemapStatusAttr();
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when target status indicates a bad mapping */
void GeneMapper::forceToUnmappedDueToTargetStatus(ResultFeatureTrees* mappedGene,
                                                  TargetStatus targetStatus) const {
    forceToUnmapped(mappedGene);
    mappedGene->setTargetStatus(targetStatus);
    mappedGene->rsetTargetStatusAttr();
}

/* check if gene contains transcripts with different mapped seqid/strand  */
bool GeneMapper::hasMixedMappedSeqStrand(const ResultFeatureTrees* mappedGene) const {
    int numTranscripts = (mappedGene->mapped != NULL) ? mappedGene->mapped->fChildren.size() : 0;
    string seqid, strand;
    for (int i = 0; i < numTranscripts; i++) {
        const FeatureNode* transcriptTree = mappedGene->mapped->fChildren[i];
        const GxfFeature* transcript = transcriptTree->fFeature;
        if (seqid == "") {
            seqid = transcript->fSeqid;
            strand = transcript->fStrand;
        } else if ((transcript->fSeqid != seqid) or (transcript->fStrand != strand)) {
            return true;
        }
    }
    return false;
}

/* check if gene contains transcripts of a gene have a target status of nonoverlap */
bool GeneMapper::hasTargetStatusNonOverlap(const ResultFeatureTrees* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return false;
    }
    if (mappedGene->mapped->fTargetStatus == TARGET_STATUS_NONOVERLAP) {
        // handles weird case of ENSG00000239810 were gene doesn't overlap
        // however the only transcript is `new'.
        return true;
    }
    for (int i = 0; i < mappedGene->mapped->fChildren.size(); i++) {
        if (mappedGene->mapped->fChildren[i]->fTargetStatus == TARGET_STATUS_NONOVERLAP) {
            return true;
        }
    }
    return false;
}

/* check if gene transcripts cause the gene size to expand beyond a
 * threshold*/
bool GeneMapper::hasExcessiveSizeChange(const ResultFeatureTrees* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return 0;  // not mapped
    } else {
        int srcGeneLength = mappedGene->src->fFeature->size();
        int mappedGeneLength = mappedGene->mapped->fFeature->size();
        return (abs(float(mappedGeneLength-srcGeneLength)) / float(srcGeneLength)) > geneExpansionThreshold;
    }
}

/* set the number of mappings of the gene  */
void GeneMapper::setNumGeneMappings(FeatureNode* mappedGeneTree) const {
    for (int i = 0; i < mappedGeneTree->fChildren.size(); i++) {
        mappedGeneTree->fNumMappings = max(mappedGeneTree->fNumMappings, mappedGeneTree->fChildren[i]->fNumMappings);
    }
}

/* should we substitute target version of gene?  */
bool GeneMapper::shouldSubstituteTarget(const ResultFeatureTrees* mappedGene) const {
    // requested and a mis-mapped gene or the right biotype to a mapped sequence
    if (not ((fSubstituteTargetVersion.size() > 0)
             and ((mappedGene->getTargetStatus() == TARGET_STATUS_NONOVERLAP)
                  or (mappedGene->getTargetStatus() == TARGET_STATUS_LOST)))) {
        return false; // not substituting or not right target status
    }
    const FeatureNode* targetGene = getTargetAnnotationNode(mappedGene->src);
    if (targetGene == NULL) {
        // this should not have happened, but it does because of a transcript id
        // ENST00000426406 incorrectly being moved to a different gene
        return false;
    }
    if (not isSrcSeqInMapping(targetGene->fFeature)) {
        return false;  // sequence not being mapped (moved chroms)
    }
    // FIXME: allow for pseudogene compatiblity
    return (targetGene->fFeature->getTypeBiotype() == mappedGene->src->fFeature->getTypeBiotype());
    return false;
}

/* copy target gene to use instead of a mapping  */
void GeneMapper::substituteTarget(ResultFeatureTrees* mappedGene) const {
    const FeatureNode* targetGene = getTargetAnnotationNode(mappedGene->src);
    mappedGene->target = targetGene->clone(mappedGene->src->fFeature->getFormat());
    // Set flag in both trees
    mappedGene->target->rsetSubstitutedMissingTargetAttr(fSubstituteTargetVersion);
    if (mappedGene->mapped != NULL) {
        mappedGene->mapped->rsetSubstitutedMissingTargetAttr(fSubstituteTargetVersion);
    }
    if (mappedGene->unmapped != NULL) {
        mappedGene->unmapped->rsetSubstitutedMissingTargetAttr(fSubstituteTargetVersion);
    }
}

/* update gene bounds given a mapped transcript */
void GeneMapper::updateMappedGeneBounds(const FeatureNode* mappedTranscript,
                                        string& seqid, string& strand,
                                        int& start, int& end) const {
    const GxfFeature* transcript = mappedTranscript->fFeature;
    assert(transcript->fType == GxfFeature::TRANSCRIPT);
    if (seqid == "") {
        // first
        seqid = transcript->fSeqid;
        strand = transcript->fStrand;
        start = transcript->fStart;
        end = transcript->fEnd;
    } else {
        start = min(start, transcript->fStart);
        end = max(end, transcript->fEnd);
    }
}

/* If there are any mapped transcripts of a gene, add a gene
 * record for the bounds */
FeatureNode* GeneMapper::buildMappedGeneFeature(const FeatureNode* srcGeneTree,
                                                ResultFeatureTreesVector& mappedTranscripts) const {
    // calculate bounds
    string seqid, strand;
    int start = 0, end = 0;
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].mapped != NULL) {
            updateMappedGeneBounds(mappedTranscripts[i].mapped, seqid, strand, start, end);
        }
    }
    FeatureNode* mappedGene = FeatureMapper::mapBounding(srcGeneTree, seqid, start-1, end, strand);  // takes zero based
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].mapped != NULL) {
            FeatureMapper::updateParent(mappedGene, mappedTranscripts[i].mapped);
        }
    }
    return mappedGene;
}

/* If there are any unmapped transcripts of a gene, add a gene
 * record for the bounds */
FeatureNode* GeneMapper::buildUnmappedGeneFeature(const FeatureNode* srcGeneTree,
                                                  ResultFeatureTreesVector& mappedTranscripts) const {
    FeatureNode* unmappedGene = FeatureMapper::mapBounding(srcGeneTree);
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].unmapped != NULL) {
            FeatureMapper::updateParent(unmappedGene, mappedTranscripts[i].unmapped);
        }
    }
    return unmappedGene;
}

/* Build gene features */
ResultFeatureTrees GeneMapper::buildGeneFeature(const FeatureNode* srcGeneTree,
                                                ResultFeatureTreesVector& mappedTranscripts) const {
    ResultFeatureTrees mappedGene(srcGeneTree);
    if (mappedTranscripts.haveMapped()) {
        mappedGene.mapped = buildMappedGeneFeature(srcGeneTree, mappedTranscripts);
        
    }
    if (mappedTranscripts.haveUnmapped()) {
        mappedGene.unmapped = buildUnmappedGeneFeature(srcGeneTree, mappedTranscripts);
    }
    if (mappedGene.mapped != NULL) {
        setNumGeneMappings(mappedGene.mapped);
    }
    TargetStatus targetStatus = getTargetAnnotationStatus(&mappedGene);
    mappedGene.setTargetStatus(targetStatus); // on gene only
    if ((targetStatus == TARGET_STATUS_LOST) && hasTargetStatusNonOverlap(&mappedGene)) {
        // lost because at least one transcript doesn't overlap the target
        mappedGene.rsetTargetStatus(TARGET_STATUS_NONOVERLAP); // gene and transcripts
    }
    return mappedGene;
}

/* output GFF3 mapped ##sequence-region if not already written */
void GeneMapper::outputMappedSeqRegionIfNeed(const FeatureNode* geneTree,
                                             ostream& mappedGxfFh) {
    if (geneTree->fFeature->getFormat() == GFF3_FORMAT) {
        const string& mappedSeqId = geneTree->fFeature->fSeqid;
        if (not checkRecordSeqRegionWritten(mappedSeqId)) {
            mappedGxfFh << "##sequence-region " << mappedSeqId << " 1 "
                        << fGenomeTransMap->getTargetSeqSize(mappedSeqId) << endl;
        }
    }
}

/*
 * recursive output of a GxF feature tree
 */
void GeneMapper::outputFeature(const FeatureNode* featureNode,
                               ostream& gxfFh) const {
    gxfFh << featureNode->fFeature->toString() << endl;
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        outputFeature(featureNode->fChildren[i], gxfFh);
    }
}

/* output genes */
void GeneMapper::outputFeatures(const ResultFeatureTrees& mappedGene,
                                ostream& mappedGxfFh,
                                ostream& unmappedGxfFh) {
    // either one of target or mapped is written
    if (mappedGene.target != NULL) {
        outputMappedSeqRegionIfNeed(mappedGene.target, mappedGxfFh);
        outputFeature(mappedGene.target, mappedGxfFh);
    } else if (mappedGene.mapped) {
        outputMappedSeqRegionIfNeed(mappedGene.mapped, mappedGxfFh);
        outputFeature(mappedGene.mapped, mappedGxfFh);
    }
    if (mappedGene.unmapped) {
        outputFeature(mappedGene.unmapped, unmappedGxfFh);
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
const GxfFeature* GeneMapper::getTargetAnnotation(const FeatureNode* featureNode) const {
    const FeatureNode* targetNode = getTargetAnnotationNode(featureNode);
    if (targetNode == NULL) {
        return NULL;
    } else {
        return targetNode->fFeature;
    }
}

/* get target annotation node for a feature, if available */
const FeatureNode* GeneMapper::getTargetAnnotationNode(const FeatureNode* featureNode) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    }
    // try id, then name
    const FeatureNode* targetNode = fTargetAnnotations->getFeatureNodeById(featureNode->fFeature->getTypeId(),
                                                                           featureNode->fFeature->fSeqid);
    if (targetNode == NULL) {
        targetNode = fTargetAnnotations->getFeatureNodeByName(featureNode->fFeature->getTypeName(),
                                                              featureNode->fFeature->fSeqid);
    }
    return targetNode;
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
TargetStatus GeneMapper::getTargetAnnotationStatus(const ResultFeatureTrees* mappedFeature) const {
    if (fTargetAnnotations == NULL) {
        return TARGET_STATUS_NA;
    }
    const GxfFeature* targetFeature = getTargetAnnotation(mappedFeature->src);
    if (targetFeature == NULL) {
        return TARGET_STATUS_NEW;
    }
    if (mappedFeature->mapped == NULL) {
        return TARGET_STATUS_LOST;
    }
    if (mappedFeature->mapped->fFeature->overlaps(targetFeature)) {
        return TARGET_STATUS_OVERLAP;
    } else {
        return TARGET_STATUS_NONOVERLAP;
    }
}

/* If target gene annotations are available, get biotype of target feature */
const string& GeneMapper::getTargetAnnotationBiotype(const ResultFeatureTrees* mappedFeature) const {
    const GxfFeature* targetFeature = getTargetAnnotation(mappedFeature->src);
    if (targetFeature == NULL) {
        return emptyString;
    } else {
        return targetFeature->getTypeBiotype();
    }
}

/* output info on one bounding feature */
void GeneMapper::outputFeatureInfo(const ResultFeatureTrees* mappedTree,
                                   bool substituteTarget,
                                   ostream& mappingInfoFh) const {
    const GxfFeature* srcFeature = mappedTree->src->fFeature;
    mappingInfoFh << srcFeature->getTypeId() << "\t"
                  << srcFeature->getTypeName() << "\t"
                  << srcFeature->fType << "\t"
                  << srcFeature->getTypeBiotype() << "\t"
                  << srcFeature->fSource << "\t"
                  << srcFeature->fSeqid << "\t"
                  << srcFeature->fStart << "\t"
                  << srcFeature->fEnd << "\t"
                  << srcFeature->fStrand << "\t";
    if ((mappedTree->mapped != NULL) && !substituteTarget) {
        const GxfFeature* mappedFeature = mappedTree->mapped->fFeature;
        mappingInfoFh << mappedFeature->fSeqid << "\t"
                      << mappedFeature->fStart << "\t"
                      << mappedFeature->fEnd << "\t"
                      << mappedFeature->fStrand << "\t";
    } else {
        mappingInfoFh << "\t0\t0\t.\t";
    }
    mappingInfoFh << remapStatusToStr(mappedTree->getRemapStatus()) << "\t"
                  << mappedTree->getNumMappings() << "\t"
                  << targetStatusToStr(mappedTree->getTargetStatus()) << "\t"
                  << getTargetAnnotationBiotype(mappedTree) << "\t"
                  << (substituteTarget ? fSubstituteTargetVersion : emptyString)
                  << endl;
}

/* find transcripts for a give source transcript and output mappings information */
void GeneMapper::outputTranscriptInfo(const ResultFeatureTrees* mappedGene,
                                      bool substituteTarget,
                                      const FeatureNode* srcTranscript,
                                      ostream& mappingInfoFh) const {
    ResultFeatureTrees transcript(srcTranscript);
    // find corresponding transcripts in each tree
    // FIXME: maybe we should keep them!!
    if (mappedGene->mapped != NULL) {
        transcript.mapped = findMatchingBoundingNode(mappedGene->mapped->fChildren, srcTranscript);
    }
    if (mappedGene->unmapped != NULL) {
        transcript.unmapped = findMatchingBoundingNode(mappedGene->unmapped->fChildren, srcTranscript);
    }
    if (mappedGene->target != NULL) {
        transcript.target = findMatchingBoundingNode(mappedGene->target->fChildren, srcTranscript);
    }
    outputFeatureInfo(&transcript, substituteTarget, mappingInfoFh);
}

/*
 * Output information about gene mapping
 */
void GeneMapper::outputInfo(const ResultFeatureTrees* mappedGene,
                            ostream& mappingInfoFh) const {
    // not all transcripts maybe not be substituted, so chec at gene level
    bool substituteTarget = mappedGene->target != NULL;
    outputFeatureInfo(mappedGene, substituteTarget, mappingInfoFh);
    // transcripts
    for (int i = 0; i < mappedGene->src->fChildren.size(); i++) {
        outputTranscriptInfo(mappedGene, substituteTarget, mappedGene->src->fChildren[i], mappingInfoFh);
    }
}

/*
 * Check for and handle problematic cases after mapping gene.
 * return true if gene is ok, false if force to unmapped.
 */
void GeneMapper::processGeneLevelMapping(ResultFeatureTrees* mappedGene) {
    if (fSkipAutomaticNonCoding and 
        isAutomaticSmallNonCodingGene(mappedGene->src)) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_AUTOMATIC_NON_CODING);
    } else if (hasMixedMappedSeqStrand(mappedGene)) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_GENE_CONFLICT);
    } else if (hasExcessiveSizeChange(mappedGene)) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_GENE_SIZE_CHANGE);
    } else if (hasTargetStatusNonOverlap(mappedGene)) {
        forceToUnmappedDueToTargetStatus(mappedGene, TARGET_STATUS_NONOVERLAP);
    }
}

/* set gene-level attributes after all mapping decisions have
 * been made */
void GeneMapper::setGeneLevelMappingAttributes(ResultFeatureTrees* mappedGene) {
    mappedGene->setBoundingFeatureRemapStatus(isSrcSeqInMapping(mappedGene->src));
    mappedGene->rsetRemapStatusAttr();
    mappedGene->setNumMappingsAttr();
    mappedGene->rsetTargetStatusAttr();
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
    FeatureNode* srcGeneTree = GeneTree::geneTreeFactory(gxfParser, geneFeature);
    ResultFeatureTreesVector mappedTranscripts = processTranscripts(srcGeneTree, transcriptPslFh);
    ResultFeatureTrees mappedGene = buildGeneFeature(srcGeneTree, mappedTranscripts);
    setGeneLevelMappingAttributes(&mappedGene);
    processGeneLevelMapping(&mappedGene);

    // must be done after forcing status above
    if (shouldSubstituteTarget(&mappedGene)) {
        substituteTarget(&mappedGene);
    }
    outputFeatures(mappedGene, mappedGxfFh, unmappedGxfFh);
    outputInfo(&mappedGene, mappingInfoFh);
    mappedGene.free();
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

/*
 * copy a target gene annotation that was skipped for mapping
 */
void GeneMapper::copySkippedTargetGene(GxfFormat gxfFormat,
                                       const FeatureNode* targetGeneNode,
                                       ostream& mappedGxfFh,
                                       ostream& mappingInfoFh) {
#if 0 // FIXME
    FeatureNode* copiedGeneNode = targetGeneNode->clone(gxfFormat);
    outputInfo(copiedGeneNode, mappingInfoFh);
#endif
}

/*
 * copy target annotations that are skipped for mapping
 */
void GeneMapper::copySkippedTargetGenes(GxfFormat gxfFormat,
                                        ostream& mappedGxfFh,
                                        ostream& mappingInfoFh) {
    const FeatureNodeVector& genes = fTargetAnnotations->getGenes();
    for (int iGene = 0; iGene < genes.size(); iGene++) {
        if (isAutomaticSmallNonCodingGene(genes[iGene])) {
            copySkippedTargetGene(gxfFormat, genes[iGene], mappedGxfFh, mappingInfoFh);
        }
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

    if (fSkipAutomaticNonCoding and (fTargetAnnotations != NULL)) {
        copySkippedTargetGenes(gxfParser->getGxfFormat(), mappedGxfFh, mappingInfoFh);
    }
}

