#include "geneMapper.hh"
#include "featureMapper.hh"
#include "jkinclude.hh"
#include "gxfRecord.hh"
#include "bedMap.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "transcriptMapper.hh"
#include "annotationSet.hh"
#include "featureTreePolish.hh"
#include "globals.hh"
#include "gxfIO.hh"


/* fraction of gene expansion that causes a rejection */
const float geneExpansionThreshold = 0.50;

/*  mapinfo TSV headers, terminated by NULL */
static const char* mappingInfoHeaders[] = {
    "id", "name", "type", "biotype", "source",
    "srcChrom", "srcStart", "srcEnd", "srcStrand",
    "mappedId", "mappedChrom", "mappedStart", "mappedEnd",
    "mappedStrand", "mappingStatus", "numMappings",
    "targetStatus", "targetBiotype", "targetSubst", NULL
};

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const Feature* feature) const {
    return fGenomeTransMap->haveQuerySeq(feature->getSeqid());
}

/* record gene or transcript as being mapped or substituted */
void GeneMapper::recordMapped(const Feature* feature) {
    fMappedIdsNames.insert(getBaseId(feature->getTypeId()));
    // N.B. gene names with `.' are not always a version
    fMappedIdsNames.insert(feature->getTypeName());
    if (feature->getHavanaTypeId() != "") {
        fMappedIdsNames.insert(getBaseId(feature->getHavanaTypeId()));
    }
}

/* record gene and it's transcript as being mapped */
void GeneMapper::recordGeneMapped(const Feature* gene) {
    recordMapped(gene);
    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        recordMapped(gene->getChild(i));
    }
}

/* check if gene or transcript have been mapped */
bool GeneMapper::checkMapped(const Feature* feature) const {
    // N.B.  Don't use gene name, as some small non-coding genes has the same
    // name for multiple instances
    if (fMappedIdsNames.find(getBaseId(feature->getTypeId())) != fMappedIdsNames.end()) {
        return true;
    }
    if (feature->getHavanaTypeId() != "") {
        if (fMappedIdsNames.find(getBaseId(feature->getHavanaTypeId())) != fMappedIdsNames.end()) {
            return true;
        }
    }
    return false;
}

/* check if gene or it's transcript have been mapped */
bool GeneMapper::checkGeneMapped(const Feature* gene) const {
    if (checkMapped(gene)) {
        return true;
    }
    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        if (checkMapped(gene->getChild(i))) {
            return true;
        }
    }
    return false;
}

/* process one transcript */
ResultFeatures GeneMapper::processTranscript(const Feature* transcript,
                                                 ostream* transcriptPslFh) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcript, fTargetAnnotations,
                                      isSrcSeqInMapping(transcript), transcriptPslFh);
    ResultFeatures mappedTranscript = transcriptMapper.mapTranscriptFeatures(transcript);
    TargetStatus targetStatus = getTargetAnnotationStatus(&mappedTranscript);
    mappedTranscript.setTargetStatus(targetStatus);
    return mappedTranscript;
}

/* process all transcripts of gene. */
ResultFeaturesVector GeneMapper::processTranscripts(const Feature* gene,
                                                    ostream* transcriptPslFh) const {
    ResultFeaturesVector mappedTranscripts;
    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        const Feature* transcript = gene->getChild(i);
        if (transcript->getType() != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcript->toString());
        }
        mappedTranscripts.push_back(processTranscript(transcript, transcriptPslFh));
    }
    return mappedTranscripts;
}

/* find a matching gene or transcript given by id */
Feature* GeneMapper::findMatchingBoundingFeature(const FeatureVector& features,
                                                 const Feature* feature) const {
    assert(feature->isGeneOrTranscript());
    for (int i = 0; i < features.size(); i++)  {
        assert(features[i]->isGeneOrTranscript());
        if (getPreMappedId(features[i]->getTypeId()) == getPreMappedId(feature->getTypeId())) {
            return features[i];  // found match!
        }
    }
    return NULL;
}

/* copy mapping metadata */
void GeneMapper::copyMappingMetadata(const Feature* origFeature,
                                     Feature* newFeature) const {
    newFeature->setRemapStatus(origFeature->getRemapStatus());
    newFeature->setTargetStatus(origFeature->getTargetStatus());
    newFeature->setNumMappings(max(origFeature->getNumMappings(), newFeature->getNumMappings()));
}

/* copy gene and transcript attributes to a fresh unmapped clone,
 * only at the transcript and gene levels. */
void GeneMapper::copyGeneMetadata(const Feature* origGene,
                                  Feature* newGene) const {
    copyMappingMetadata(origGene, newGene);

    // copy transcript metadata
    for (size_t i = 0; i < newGene->getChildren().size(); i++) {
        Feature* newTranscript = newGene->getChild(i);
        const Feature* origTranscript = findMatchingBoundingFeature(origGene->getChildren(), newTranscript);
        if (origTranscript != NULL) {
            copyMappingMetadata(origTranscript, newTranscript);
        }
    }
}

/* force to unmapped */
void GeneMapper::forceToUnmapped(ResultFeatures* mappedGene) const {
    // create copy and update attributes at gene/transcript level before freeing.
    Feature* srcCopy = mappedGene->src->cloneTree();
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
void GeneMapper::forceToUnmappedDueToRemapStatus(ResultFeatures* mappedGene,
                                                 RemapStatus remapStatus) const {
    if (gVerbose) {
        cerr << "forceToUnmappedDueToRemapStatus " << mappedGene->src->getTypeId()
             << " " << remapStatusToStr(remapStatus) << endl;
    }
    forceToUnmapped(mappedGene);
    mappedGene->rsetRemapStatus(remapStatus);
    mappedGene->rsetRemapStatusAttr();
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when target status indicates a bad mapping */
void GeneMapper::forceToUnmappedDueToTargetStatus(ResultFeatures* mappedGene,
                                                  TargetStatus targetStatus) const {
    if (gVerbose) {
        cerr << "forceToUnmappedDueToTargetStatus " << mappedGene->src->getTypeId()
             << " " << targetStatusToStr(targetStatus) << endl;
    }
    forceToUnmapped(mappedGene);
    mappedGene->setTargetStatus(targetStatus);
    mappedGene->rsetTargetStatusAttr();
}

/* check if gene contains transcripts with different mapped seqid/strand  */
bool GeneMapper::hasMixedMappedSeqStrand(const ResultFeatures* mappedGene) const {
    int numTranscripts = (mappedGene->mapped != NULL) ? mappedGene->mapped->getChildren().size() : 0;
    string seqid, strand;
    for (int i = 0; i < numTranscripts; i++) {
        const Feature* transcript = mappedGene->mapped->getChild(i);
        if (seqid == "") {
            seqid = transcript->getSeqid();
            strand = transcript->getStrand();
        } else if ((transcript->getSeqid() != seqid) or (transcript->getStrand() != strand)) {
            return true;
        }
    }
    return false;
}

/* check if gene contains transcripts of a gene have a target status of nonoverlap */
bool GeneMapper::hasTargetStatusNonOverlap(const ResultFeatures* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return false;
    }
    if (mappedGene->mapped->getTargetStatus() == TARGET_STATUS_NONOVERLAP) {
        // handles weird case of ENSG00000239810 were gene doesn't overlap
        // however the only transcript is `new'.
        return true;
    }
    for (int i = 0; i < mappedGene->mapped->getChildren().size(); i++) {
        if (mappedGene->mapped->getChild(i)->getTargetStatus() == TARGET_STATUS_NONOVERLAP) {
            return true;
        }
    }
    return false;
}

/* check if gene transcripts cause the gene size to expand beyond a
 * threshold*/
bool GeneMapper::hasExcessiveSizeChange(const ResultFeatures* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return 0;  // not mapped
    } else {
        int srcGeneLength = mappedGene->src->size();
        int mappedGeneLength = mappedGene->mapped->size();
        return (abs(float(mappedGeneLength-srcGeneLength)) / float(srcGeneLength)) > geneExpansionThreshold;
    }
}

/* set the number of mappings of the gene  */
void GeneMapper::setNumGeneMappings(Feature* mappedGeneTree) const {
    for (int i = 0; i < mappedGeneTree->getChildren().size(); i++) {
        mappedGeneTree->setNumMappings(max(mappedGeneTree->getNumMappings(), mappedGeneTree->getChild(i)->getNumMappings()));
    }
}

/*
 * Check for pathological case: gene names have been moved to different ids
 * and mapped gene ids doesn't exist, so match was made on gene name.  This
 * result in a gene_id being include twice. Tested for by CTC-559E9.12 and
 * BNIP3P9 entries.  if target gene_id doesn't match source gene_id and target
 * gene_id is in source mappings, we should not substitute the target.
 */
bool GeneMapper::checkForPathologicalGeneRename(const ResultFeatures* mappedGene,
                                                const Feature* targetGene) const {
    return (getBaseId(mappedGene->src->getTypeId()) != getBaseId(targetGene->getTypeId()))
        and (fSrcAnnotations->getFeatureById(targetGene->getTypeId(), targetGene->getSeqid()) != NULL);
}

/* should we substitute target version of gene?  */
bool GeneMapper::shouldSubstituteTarget(const ResultFeatures* mappedGene) const {
    // mis-mapped gene or the right biotype to a mapped sequence
    if (fSubstituteTargetVersion.size() == 0) {
        return false; // not substituting
    }
    const Feature* targetGene = getTargetAnnotation(mappedGene->src);
    if (targetGene == NULL) {
        // this should not have happened, but it does because of a transcript id
        // ENST00000426406 incorrectly being moved to a different gene
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << mappedGene->src->getTypeId()
                 << " missing target gene due to gene id rename" << endl;
        }
        return false;
    }
    if (not isSrcSeqInMapping(targetGene)) {
        return false;  // sequence not being mapped (moved chroms)
    }

    if (checkForPathologicalGeneRename(mappedGene, targetGene)) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << mappedGene->src->getTypeId()
                 << "  " << mappedGene->src->getTypeName()
                 << " pathological gene rename" << endl;
        }
        return false;
    }
    
    // allow for pseudogene biotype compatiblity between less and more
    // specific pseudogene biotypes
    if ((targetGene->getTypeBiotype() == mappedGene->src->getTypeBiotype())
        or (targetGene->isPseudogene() == mappedGene->src->isPseudogene())) {
        // must check to make sure it wasn't already mapped due to gene id/name pairing incorrectly changing
        bool alreadyMapped = checkGeneMapped(targetGene);
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: " << !alreadyMapped << ": mapped: " << mappedGene->src->getTypeId()
                 << " substitute target: " << targetGene->getTypeId() << " alreadyMapped: " << alreadyMapped << endl;
        }
        return !alreadyMapped;
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

/* copy target gene to use instead of a mapping  */
void GeneMapper::substituteTarget(ResultFeatures* mappedGene) {
    mappedGene->target = getTargetAnnotation(mappedGene->src)->cloneTree();

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
void GeneMapper::updateMappedGeneBounds(const Feature* mappedTranscript,
                                        string& seqid, string& strand,
                                        int& start, int& end) const {
    assert(mappedTranscript->getType() == GxfFeature::TRANSCRIPT);
    if (seqid == "") {
        // first
        seqid = mappedTranscript->getSeqid();
        strand = mappedTranscript->getStrand();
        start = mappedTranscript->getStart();
        end = mappedTranscript->getEnd();
    } else {
        start = min(start, mappedTranscript->getStart());
        end = max(end, mappedTranscript->getEnd());
    }
}

/* If there are any mapped transcripts of a gene, add a gene
 * record for the bounds */
Feature* GeneMapper::buildMappedGeneFeature(const Feature* srcGeneTree,
                                                ResultFeaturesVector& mappedTranscripts) const {
    // calculate bounds
    string seqid, strand;
    int start = 0, end = 0;
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].mapped != NULL) {
            updateMappedGeneBounds(mappedTranscripts[i].mapped, seqid, strand, start, end);
        }
    }
    Feature* mappedGene = FeatureMapper::mapBounding(srcGeneTree, seqid, start-1, end, strand);  // takes zero based
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].mapped != NULL) {
            FeatureMapper::updateParent(mappedGene, mappedTranscripts[i].mapped);
        }
    }
    return mappedGene;
}

/* If there are any unmapped transcripts of a gene, add a gene
 * record for the bounds */
Feature* GeneMapper::buildUnmappedGeneFeature(const Feature* srcGeneTree,
                                                  ResultFeaturesVector& mappedTranscripts) const {
    Feature* unmappedGene = FeatureMapper::mapBounding(srcGeneTree);
    for (int i = 0; i < mappedTranscripts.size(); i++) {
        if (mappedTranscripts[i].unmapped != NULL) {
            FeatureMapper::updateParent(unmappedGene, mappedTranscripts[i].unmapped);
        }
    }
    return unmappedGene;
}

/* Build gene features */
ResultFeatures GeneMapper::buildGeneFeature(const Feature* srcGeneTree,
                                                ResultFeaturesVector& mappedTranscripts) const {
    ResultFeatures mappedGene(srcGeneTree);
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
void GeneMapper::saveMapped(ResultFeatures& mappedGene,
                            AnnotationSet& mappedSet) {
    // either one of target or mapped is saved
    if (mappedGene.target != NULL) {
        recordGeneMapped(mappedGene.target);
        mappedSet.addGene(mappedGene.target);
        mappedGene.target = NULL;
    } else if (mappedGene.mapped != NULL) {
        recordGeneMapped(mappedGene.mapped);
        mappedSet.addGene(mappedGene.mapped);
        mappedGene.mapped = NULL;
    }
}

/* save unmapped gene features  */
void GeneMapper::saveUnmapped(ResultFeatures& mappedGene,
                              AnnotationSet& unmappedSet) {
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

/* get target annotation  for a feature, if available */
const Feature* GeneMapper::getTargetAnnotation(const Feature* feature) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    }
    // try id, havana id, then name
    const Feature* targetFeature = fTargetAnnotations->getFeatureById(getPreMappedId(feature->getTypeId()),
                                                                      feature->getSeqid());
    if ((targetFeature == NULL) and (feature->getHavanaTypeId() != "")) {
        targetFeature = fTargetAnnotations->getFeatureByName(getPreMappedId(feature->getHavanaTypeId()),
                                                             feature->getSeqid());
    }
    if (targetFeature == NULL) {
        targetFeature = fTargetAnnotations->getFeatureByName(getPreMappedId(feature->getTypeName()),
                                                             feature->getSeqid());
    }
    return targetFeature;
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
TargetStatus GeneMapper::getTargetAnnotationStatus(const ResultFeatures* mappedFeature) const {
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
    if (mappedFeature->mapped->overlaps(targetFeature)) {
        return TARGET_STATUS_OVERLAP;
    } else {
        return TARGET_STATUS_NONOVERLAP;
    }
}

/* If target gene annotations are available, get biotype of target feature */
const string& GeneMapper::getTargetAnnotationBiotype(const ResultFeatures* mappedFeature) const {
    const GxfFeature* targetFeature = getTargetAnnotation(mappedFeature->src);
    if (targetFeature == NULL) {
        return emptyString;
    } else {
        return targetFeature->getTypeBiotype();
    }
}

/* output info on one bounding feature */
void GeneMapper::outputFeatureInfo(const ResultFeatures* mappedTree,
                                   bool substituteTarget,
                                   ostream& mappingInfoFh) const {
    const GxfFeature* srcFeature = mappedTree->src;
    mappingInfoFh << srcFeature->getTypeId() << "\t"
                  << srcFeature->getTypeName() << "\t"
                  << srcFeature->getType() << "\t"
                  << srcFeature->getTypeBiotype() << "\t"
                  << srcFeature->getSource() << "\t"
                  << srcFeature->getSeqid() << "\t"
                  << srcFeature->getStart() << "\t"
                  << srcFeature->getEnd() << "\t"
                  << srcFeature->getStrand() << "\t";
    if ((mappedTree->mapped != NULL) && !substituteTarget) {
        const GxfFeature* mappedFeature = mappedTree->mapped;
        mappingInfoFh << mappedFeature->getTypeId() << "\t"
                      << mappedFeature->getSeqid() << "\t"
                      << mappedFeature->getStart() << "\t"
                      << mappedFeature->getEnd() << "\t"
                      << mappedFeature->getStrand() << "\t";
    } else {
        mappingInfoFh << "\t\t0\t0\t.\t";
    }
    mappingInfoFh << remapStatusToStr(mappedTree->getRemapStatus()) << "\t"
                  << mappedTree->getNumMappings() << "\t"
                  << targetStatusToStr(mappedTree->getTargetStatus()) << "\t"
                  << getTargetAnnotationBiotype(mappedTree) << "\t"
                  << (substituteTarget ? fSubstituteTargetVersion : emptyString)
                  << endl;
}

/* find transcripts for a give source transcript and output mappings information */
void GeneMapper::outputTranscriptInfo(const ResultFeatures* mappedGene,
                                      bool substituteTarget,
                                      const Feature* srcTranscript,
                                      ostream& mappingInfoFh) const {
    ResultFeatures transcript(srcTranscript);
    // find corresponding transcripts in each tree
    if (mappedGene->mapped != NULL) {
        transcript.mapped = findMatchingBoundingFeature(mappedGene->mapped->getChildren(), srcTranscript);
    }
    if (mappedGene->unmapped != NULL) {
        transcript.unmapped = findMatchingBoundingFeature(mappedGene->unmapped->getChildren(), srcTranscript);
    }
    if (mappedGene->target != NULL) {
        transcript.target = findMatchingBoundingFeature(mappedGene->target->getChildren(), srcTranscript);
    }
    outputFeatureInfo(&transcript, substituteTarget, mappingInfoFh);
}

/*
 * Output information about gene mapping
 */
void GeneMapper::outputInfo(const ResultFeatures* mappedGene,
                            ostream& mappingInfoFh) const {
    // not all transcripts maybe not be substituted, so chec at gene level
    bool substituteTarget = mappedGene->target != NULL;
    outputFeatureInfo(mappedGene, substituteTarget, mappingInfoFh);
    // transcripts
    for (int i = 0; i < mappedGene->src->getChildren().size(); i++) {
        outputTranscriptInfo(mappedGene, substituteTarget, mappedGene->src->getChild(i), mappingInfoFh);
    }
}

/*
 * Check for and handle problematic cases after mapping gene.
 * return true if gene is ok, false if force to unmapped.
 */
void GeneMapper::processGeneLevelMapping(ResultFeatures* mappedGene) {
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
void GeneMapper::setGeneLevelMappingAttributes(ResultFeatures* mappedGene) {
    mappedGene->setBoundingFeatureRemapStatus(isSrcSeqInMapping(mappedGene->src));
    mappedGene->rsetRemapStatusAttr();
    mappedGene->setNumMappingsAttr();
    mappedGene->rsetTargetStatusAttr();
}

/*
 * map and output one gene's annotations
 */
void GeneMapper::mapGene(const Feature* srcGeneTree,
                         AnnotationSet& mappedSet,
                         AnnotationSet& unmappedSet,
                         FeatureTreePolish& featureTreePolish,
                         ostream& mappingInfoFh,
                         ostream* transcriptPslFh) {
    ResultFeaturesVector mappedTranscripts = processTranscripts(srcGeneTree, transcriptPslFh);
    ResultFeatures mappedGene = buildGeneFeature(srcGeneTree, mappedTranscripts);
    setGeneLevelMappingAttributes(&mappedGene);
    processGeneLevelMapping(&mappedGene);

    // must be done after forcing status above
    if ((mappedGene.mapped == NULL) and shouldSubstituteTarget(&mappedGene)) {
        substituteTarget(&mappedGene);
    }
    if (mappedGene.mapped != NULL) {
        featureTreePolish.polishGene(mappedGene.mapped);
    }
    outputInfo(&mappedGene, mappingInfoFh);  // MUST do before saveGene, as it moved to output sets
    saveMapped(mappedGene, mappedSet);
    saveUnmapped(mappedGene, unmappedSet);
    mappedGene.free();
}

/* determine if this is a gene type that should not be mapped, returning
 * the remap status */
RemapStatus GeneMapper::getNoMapRemapStatus(const Feature* gene) const {
    if ((fUseTargetFlags & useTargetForAutoNonCoding) && gene->isAutomaticSmallNonCodingGene()) {
        return REMAP_STATUS_AUTO_SMALL_NCRNA;
    } else if ((fUseTargetFlags & useTargetForAutoGenes)  && gene->isAutomatic()) {
        return REMAP_STATUS_AUTOMATIC_GENE;
    } else if ((fUseTargetFlags & useTargetForPseudoGenes) && gene->isPseudogene()) {
        return REMAP_STATUS_PSEUDOGENE;
    } else {
        return REMAP_STATUS_NONE;
    }
} 

/*
 * check if a gene type should be mapped or targets of this type substituted.
 */
bool GeneMapper::shouldMapGeneType(const Feature* gene) const {
    return getNoMapRemapStatus(gene) == REMAP_STATUS_NONE;
}

/* is target gene in patch region? */
bool GeneMapper::inTargetPatchRegion(const Feature* targetGene) {
    return fTargetPatchMap->anyOverlap(targetGene->getSeqid(),
                                       targetGene->getStart(),
                                       targetGene->getEnd());
}

/* check to see if the target overlaps a mapped gene with sufficient similarity to
 * be considered the same annotation..  */
bool GeneMapper::checkTargetOverlappingMapped(const Feature* targetGene,
                                              AnnotationSet& mappedSet) {
    static const float minSimilarity = 0.5;
    FeatureVector overlapping = mappedSet.findOverlappingGenes(targetGene, minSimilarity,
                                                               fOnlyManualForTargetSubstituteOverlap);
    return overlapping.size() > 0;
}

/*
 * Check if a target gene should be copied.
 */
bool GeneMapper::shouldIncludeTargetGene(const Feature* targetGene,
                                         AnnotationSet& mappedSet)  {
    if (gVerbose) {
        cerr << "shouldIncludeTargetGene: " << targetGene->getTypeId()
             << "  " << targetGene->getTypeBiotype()
             << " noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(targetGene))
             << " shouldMapGeneType: " << shouldMapGeneType(targetGene)
             << " already mapped: " << checkGeneMapped(targetGene)
             << endl;
    }
    if (not shouldMapGeneType(targetGene)) {
        // biotypes not excluding from mapped, checkGeneMapped handles
        // case where biotype has changed
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: isSrcSeqInMapping:" << isSrcSeqInMapping(targetGene) << endl;
        }
        return isSrcSeqInMapping(targetGene) and not checkGeneMapped(targetGene);
    }
    if ((fUseTargetFlags & useTargetForPatchRegions) && inTargetPatchRegion(targetGene)) {
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: in patched region: "
                 << " overlaps mapping: " << checkTargetOverlappingMapped(targetGene, mappedSet) << endl;
        }
        // don't use if there is a mapped with significant overlap
        return (not checkGeneMapped(targetGene))
            and (not checkTargetOverlappingMapped(targetGene, mappedSet));
    }
    return false;
}

/*
 * copy a target gene annotation that was skipped for mapping
 */
void GeneMapper::copyTargetGene(const Feature* targetGene,
                                AnnotationSet& mappedSet,
                                ostream& mappingInfoFh) {
    ResultFeatures mappedGene;
    mappedGene.src = targetGene;
    mappedGene.target = targetGene->cloneTree();
    mappedGene.target->rsetRemapStatus(getNoMapRemapStatus(targetGene));
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
    const FeatureVector& genes = fTargetAnnotations->getGenes();
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
    FeatureTreePolish featureTreePolish(fPreviousMappedAnotations);
    
    const FeatureVector& srcGenes = fSrcAnnotations->getGenes();
    outputInfoHeader(mappingInfoFh);
    for (int i = 0; i < srcGenes.size(); i++) {
        if (gVerbose) {
            cerr << endl << "mapGxf: " << srcGenes[i]->getTypeId() << " " << srcGenes[i]->getTypeName()
                 << " shouldMapGeneType: " << shouldMapGeneType(srcGenes[i])
                 << " noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(srcGenes[i]))
                 << " " << srcGenes[i]->getTypeId() << " " << srcGenes[i]->getSource()
                 << endl;
        }
        if (shouldMapGeneType(srcGenes[i])) {
            mapGene(srcGenes[i], mappedSet, unmappedSet, featureTreePolish, mappingInfoFh, transcriptPslFh);
        }
    }
    if ((fUseTargetFlags != 0) and (fTargetAnnotations != NULL)) {
        copyTargetGenes(mappedSet, mappingInfoFh);
    }
    mappedSet.write(mappedGxfFh);
    unmappedSet.write(unmappedGxfFh);
}

