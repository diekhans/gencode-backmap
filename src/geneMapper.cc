#include "geneMapper.hh"
#include "featureMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "bedMap.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "transcriptMapper.hh"
#include "annotationSet.hh"
#include "globals.hh"


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

 /* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const GxfFeature* feature) const {
    return fGenomeTransMap->haveQuerySeq(feature->fSeqid);
}

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const FeatureNode* featureNode) const {
    return isSrcSeqInMapping(featureNode->fFeature);
}

/* record gene or transcript as being mapped */
void GeneMapper::recordMapped(const FeatureNode* featureNode) {
    fMappedIdsNames.insert(getBaseId(featureNode->getTypeId()));
    // N.B. gene names with `.' are not always a version
    fMappedIdsNames.insert(featureNode->getTypeName());
    if (featureNode->getHavanaTypeId() != "") {
        fMappedIdsNames.insert(getBaseId(featureNode->getHavanaTypeId()));
    }
}

/* record gene and it's transcript as being mapped */
void GeneMapper::recordGeneMapped(const FeatureNode* geneTree) {
    recordMapped(geneTree);
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        recordMapped(geneTree->fChildren[i]);
    }
}

/* check if gene or transcript have been mapped */
bool GeneMapper::checkMapped(const FeatureNode* featureNode) const {
    if (fMappedIdsNames.find(getBaseId(featureNode->getTypeId())) != fMappedIdsNames.end()) {
        return true;
    }
    if (fMappedIdsNames.find(featureNode->getTypeName()) != fMappedIdsNames.end()) {
        return true;
    }
    if (featureNode->getHavanaTypeId() != "") {
        if (fMappedIdsNames.find(getBaseId(featureNode->getHavanaTypeId())) != fMappedIdsNames.end()) {
            return true;
        }
    }
    return false;
}

/* check if gene or it's transcript have been mapped */
bool GeneMapper::checkGeneMapped(const FeatureNode* geneTree) const {
    if (checkMapped(geneTree)) {
        return true;
    }
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        if (checkMapped(geneTree->fChildren[i])) {
            return true;
        }
    }
    return false;
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
        if (features[i]->getTypeId() == feature->getTypeId()) {
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
    FeatureNode* srcCopy = mappedGene->src->clone();
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
    if (fSubstituteTargetVersion.size() == 0) {
        return false; // not substituting
    }
    if (not ((mappedGene->getTargetStatus() == TARGET_STATUS_NONOVERLAP)
             or (mappedGene->getTargetStatus() == TARGET_STATUS_LOST))) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << mappedGene->src->getTypeId()
                 << " wrong target status: " << targetStatusToStr(mappedGene->getTargetStatus())
                 << endl;
        }
        return false; // not right target status
    }
    const FeatureNode* targetGene = getTargetAnnotationNode(mappedGene->src);
    if (targetGene == NULL) {
        // this should not have happened, but it does because of a transcript id
        // ENST00000426406 incorrectly being moved to a different gene
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << mappedGene->src->getTypeId()
                 << " missing target gene due to gene id rename" << endl;
        }
        return false;
    }
    if (not isSrcSeqInMapping(targetGene->fFeature)) {
        return false;  // sequence not being mapped (moved chroms)
    }
    // allow for pseudogene biotype compatiblity between less and more
    // specific pseudogene biotypes
    if ((targetGene->getTypeBiotype() == mappedGene->src->getTypeBiotype())
        or (targetGene->isPseudogene() == mappedGene->src->isPseudogene())) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: true: mapped: " << mappedGene->src->getTypeId()
                 << " substitute target: " << targetGene->getTypeId() << endl;
        }
        return true;
    } else {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << mappedGene->src->getTypeId()
                 << " gene biotype change: target: " << targetGene->getTypeBiotype()
                 << " mapped: " << mappedGene->src->getTypeBiotype()
                 << endl;
        }
        return false;
    }
}

/* clone a target gene and remove 'transcript_*'' attributes that an old bug
 * put on gene records.
 */
FeatureNode* GeneMapper::cloneTargetGene(const FeatureNode* srcGene) const {
    FeatureNode* targetGene = getTargetAnnotationNode(srcGene)->clone();
    AttrVals& attrVals = targetGene->fFeature->getAttrs();
    attrVals.remove("transcript_id");
    attrVals.remove("transcript_type");
    attrVals.remove("transcript_name");
    attrVals.remove("transcript_status");
    return targetGene;
}

/* copy target gene to use instead of a mapping  */
void GeneMapper::substituteTarget(ResultFeatureTrees* mappedGene) const {
    mappedGene->target = cloneTargetGene(mappedGene->src);

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
    if (gVerbose) {
        cerr << "buildGeneFeature: "  << srcGeneTree->getTypeId()
             << " " << remapStatusToStr(mappedGene.getRemapStatus())
             << " " << targetStatusToStr(mappedGene.getTargetStatus()) << endl;
    }
    return mappedGene;
}

/* save mapped gene features  */
void GeneMapper::saveMapped(ResultFeatureTrees& mappedGene,
                            AnnotationSet& mappedSet) const {
    // either one of target or mapped is saved
    if (mappedGene.target != NULL) {
        mappedSet.addGene(mappedGene.target);
        mappedGene.target = NULL;
    } else if (mappedGene.mapped != NULL) {
        mappedSet.addGene(mappedGene.mapped);
        mappedGene.mapped = NULL;
    }
}

/* save unmapped gene features  */
void GeneMapper::saveUnmapped(ResultFeatureTrees& mappedGene,
                              AnnotationSet& unmappedSet) const {
    if (mappedGene.unmapped) {
        unmappedSet.addGene(mappedGene.unmapped);
        mappedGene.unmapped = NULL;
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
    // try id, havana id, then name
    const FeatureNode* targetNode = fTargetAnnotations->getFeatureNodeById(featureNode->getTypeId(),
                                                                           featureNode->fFeature->fSeqid);
    if ((targetNode == NULL) and (featureNode->getHavanaTypeId() != "")) {
        targetNode = fTargetAnnotations->getFeatureNodeByName(featureNode->getHavanaTypeId(),
                                                              featureNode->fFeature->fSeqid);
    }
    if (targetNode == NULL) {
        targetNode = fTargetAnnotations->getFeatureNodeByName(featureNode->getTypeName(),
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
    if (hasMixedMappedSeqStrand(mappedGene)) {
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
void GeneMapper::mapGene(const FeatureNode* srcGeneTree,
                         AnnotationSet& mappedSet,
                         AnnotationSet& unmappedSet,
                         ostream& mappingInfoFh,
                         ostream* transcriptPslFh) {
    if (gVerbose) {
        cerr << endl << "mapGene: " << srcGeneTree->getTypeId() << endl;
    }
    ResultFeatureTreesVector mappedTranscripts = processTranscripts(srcGeneTree, transcriptPslFh);
    ResultFeatureTrees mappedGene = buildGeneFeature(srcGeneTree, mappedTranscripts);
    setGeneLevelMappingAttributes(&mappedGene);
    processGeneLevelMapping(&mappedGene);

    // must be done after forcing status above
    if (shouldSubstituteTarget(&mappedGene)) {
        substituteTarget(&mappedGene);
    }
    outputInfo(&mappedGene, mappingInfoFh);  // MUST do before saveGene, as it moved to output sets
    saveMapped(mappedGene, mappedSet);
    saveUnmapped(mappedGene, unmappedSet);
    if (mappedGene.mapped != NULL) {
        recordGeneMapped(mappedGene.src);
    }
    mappedGene.free();
}

/* determine if this is a gene type that should not be mapped, returning
 * the remap status */
RemapStatus GeneMapper::getNoMapRemapStatus(const FeatureNode* geneTree) const {
    if ((fUseTargetFlags & useTargetForAutoNonCoding) && geneTree->isAutomaticSmallNonCodingGene()) {
        return REMAP_STATUS_AUTO_SMALL_NCRNA;
    } else if ((fUseTargetFlags & useTargetForAutoGenes)  && geneTree->isAutomatic()) {
        return REMAP_STATUS_AUTOMATIC_GENE;
    } else if ((fUseTargetFlags & useTargetForPseudoGenes) && geneTree->isPseudogene()) {
        return REMAP_STATUS_PSEUDOGENE;
    } else {
        return REMAP_STATUS_NONE;
    }
} 

/*
 * check if a gene type should be mapped or targets of this type substituted.
 */
bool GeneMapper::shouldMapGeneType(const FeatureNode* geneTree) const {
    return getNoMapRemapStatus(geneTree) == REMAP_STATUS_NONE;
}

/* is target gene in patch region? */
bool GeneMapper::inTargetPatchRegion(const FeatureNode* targetGene) {
    return fTargetPatchMap->anyOverlap(targetGene->fFeature->fSeqid,
                                       targetGene->fFeature->fStart,
                                       targetGene->fFeature->fEnd);
}

/* check to see if the target overlaps a mapped gene with sufficient similarity to
 * be considered the same annotation..  */
bool GeneMapper::checkTargetOverlappingMapped(const FeatureNode* targetGene,
                                              AnnotationSet& mappedSet) {
    static const float minSimilarity = 0.5;
    FeatureNodeVector overlapping = mappedSet.findOverlappingGenes(targetGene, minSimilarity);
    return overlapping.size() > 0;
}

/*
 * Check if a target gene should be copied.
 */
bool GeneMapper::shouldIncludeTargetGene(const FeatureNode* targetGene,
                                         AnnotationSet& mappedSet)  {
        if (gVerbose) {
            cerr << "shouldIncludeTargetGene: " << targetGene->getTypeId()
                 << "  " << targetGene->getTypeBiotype() << endl;
        }
    bool shouldInclude = false;
    if (not shouldMapGeneType(targetGene)) {
        // biotypes not excluding from mapped, checkGeneMapped handles
        // case where biotype has changed
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: isSrcSeqInMapping:" << isSrcSeqInMapping(targetGene)
                 << " already mapped: " << checkGeneMapped(targetGene) << endl;
        }
        return isSrcSeqInMapping(targetGene) and not checkGeneMapped(targetGene);
    }
    if ((fUseTargetFlags & useTargetForPatchRegions) && inTargetPatchRegion(targetGene)) {
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: in patched region: already mapped: " << checkGeneMapped(targetGene)
                 << " overlaps mapping: " << checkTargetOverlappingMapped(targetGene, mappedSet) << endl;
        }
        // don't use if there is a mapped with significant overlap
        return (not checkGeneMapped(targetGene))
            and (not checkTargetOverlappingMapped(targetGene, mappedSet));
    }
    
    if (gVerbose) {
        cerr << "    shouldIncludeTargetGene: " << shouldInclude  << endl;
    }
    return shouldInclude;
}

/*
 * copy a target gene annotation that was skipped for mapping
 */
void GeneMapper::copyTargetGene(const FeatureNode* targetGeneNode,
                                AnnotationSet& mappedSet,
                                ostream& mappingInfoFh) {
    ResultFeatureTrees mappedGene;
    mappedGene.src = targetGeneNode;
    mappedGene.target = targetGeneNode->clone();
    mappedGene.target->rsetRemapStatus(getNoMapRemapStatus(targetGeneNode));
    mappedGene.target->rsetRemapStatusAttr();
    mappedGene.target->rsetTargetStatusAttr();
    mappedGene.target->rsetSubstitutedMissingTargetAttr(fSubstituteTargetVersion);
    outputInfo(&mappedGene, mappingInfoFh); // MUST do before copying
    saveMapped(mappedGene, mappedSet);
    mappedGene.src = NULL; // don't free!!
}

/*
 * copy target annotations that are skipped for mapping
 */
void GeneMapper::copyTargetGenes(AnnotationSet& mappedSet,
                                 ostream& mappingInfoFh) {
    const FeatureNodeVector& genes = fTargetAnnotations->getGenes();
    for (int iGene = 0; iGene < genes.size(); iGene++) {
        if (shouldIncludeTargetGene(genes[iGene], mappedSet)) {
            copyTargetGene(genes[iGene], mappedSet, mappingInfoFh);
        }
    }
}

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfWriter& mappedGxfFh,
                        GxfWriter& unmappedGxfFh,
                        ostream& mappingInfoFh,
                        ostream* transcriptPslFh) {
    AnnotationSet mappedSet(&fGenomeTransMap->fTargetSizes);
    AnnotationSet unmappedSet(&fGenomeTransMap->fQuerySizes);
    
    const FeatureNodeVector& srcGenes = fSrcAnnotations->getGenes();
    outputInfoHeader(mappingInfoFh);
    for (int i = 0; i < srcGenes.size(); i++) {
        if (shouldMapGeneType(srcGenes[i])) {
            mapGene(srcGenes[i], mappedSet, unmappedSet, mappingInfoFh, transcriptPslFh);
        }
    }
    if ((fUseTargetFlags != 0) and (fTargetAnnotations != NULL)) {
        copyTargetGenes(mappedSet, mappingInfoFh);
    }
    mappedSet.write(mappedGxfFh);
    unmappedSet.write(unmappedGxfFh);
}

