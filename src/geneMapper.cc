#include "geneMapper.hh"
#include "featureMapper.hh"
#include "jkinclude.hh"
#include "feature.hh"
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
    "geneNum", "recType", "featType",
    "featId", "featOttId", "featName", "featBiotype", "featChrom", "featStart", "featEnd", "featStrand",
    "mappingStatus", "mappingCount", "targetStatus", NULL
};

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

/* output info record */
void GeneMapper::outputInfo(const string& recType,
                            const string& featType,
                            const Feature* feature,
                            RemapStatus mappingStatus,
                            int mappingCount,
                            TargetStatus targetStatus,
                            ostream& mappingInfoFh) const {
    mappingInfoFh << fCurrentGeneNum << "\t"
                  << recType << "\t"
                  << featType << "\t"
                  << feature->getTypeId() << "\t"
                  << feature->getHavanaTypeId() << "\t"
                  << feature->getTypeName() << "\t"
                  << feature->getTypeBiotype() << "\t"
                  << feature->getSeqid() << "\t"
                  << feature->getStart() << "\t"
                  << feature->getEnd() << "\t"
                  << feature->getStrand() << "\t"
                  << remapStatusToStr(mappingStatus) << "\t"
                  << mappingCount << "\t"
                  << targetStatusToStr(targetStatus)
                  << endl;
}

/*
 * Output information about source gene that was  mapped or failed mapping.
 */
void GeneMapper::outputSrcGeneInfo(const ResultFeatures* mappedGene,
                                   ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    const Feature* srcGene = mappedGene->src;
    outputInfo("mapSrc", "gene", srcGene, srcGene->getRemapStatus(), 0, srcGene->getTargetStatus(), mappingInfoFh);
    for (int i = 0; i < srcGene->getChildren().size(); i++) {
        const Feature* srcTrans = srcGene->getChild(i);
        outputInfo("mapSrc", "trans", srcTrans, srcTrans->getRemapStatus(), 0, srcTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene that was mapped or partially mapped.
 */
void GeneMapper::outputMappedGeneInfo(const ResultFeatures* mappedGene,
                                      ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    // not all transcripts maybe not be mapped
    const Feature* mapGene = mappedGene->mapped;
    outputInfo("map", "gene", mapGene, mapGene->getRemapStatus(), mapGene->getNumMappings(), mapGene->getTargetStatus(), mappingInfoFh);

    // transcripts
    for (int i = 0; i < mapGene->getChildren().size(); i++) {
        const Feature* mapTrans = mapGene->getChild(i);
        outputInfo("map", "trans", mapTrans, mapTrans->getRemapStatus(), mapGene->getNumMappings(), mapTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene that was unmapped or partially unmapped.
 */
void GeneMapper::outputUnmappedGeneInfo(const ResultFeatures* mappedGene,
                                        ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    // not all transcripts maybe not be mapped
    const Feature* unmapGene = mappedGene->unmapped;
    outputInfo("unmap", "gene", unmapGene, unmapGene->getRemapStatus(), unmapGene->getNumMappings(), unmapGene->getTargetStatus(), mappingInfoFh);

    // transcripts
    for (int i = 0; i < unmapGene->getChildren().size(); i++) {
        const Feature* unmapTrans = unmapGene->getChild(i);
        outputInfo("unmap", "trans", unmapTrans, unmapTrans->getRemapStatus(), unmapTrans->getNumMappings(), unmapTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene where target is used by default or substituted
 */
void GeneMapper::outputTargetGeneInfo(const ResultFeatures* mappedGene,
                                      const string& targetAction,
                                      ostream& mappingInfoFh) const {
    const Feature* targetGene = mappedGene->target;
    outputInfo(targetAction, "gene", targetGene, targetGene->getRemapStatus(), 0, targetGene->getTargetStatus(), mappingInfoFh);
    for (int i = 0; i < targetGene->getChildren().size(); i++) {
        const Feature* targetTrans = targetGene->getChild(i);
        outputInfo(targetAction, "trans", targetTrans, targetTrans->getRemapStatus(), 0, targetTrans->getTargetStatus(), mappingInfoFh);
    }
}

/* get string describing gene or other feature */
string GeneMapper::featureDesc(const Feature* gene) const {
    return "(" + gene->getTypeId() + " " + gene->getTypeName() + " " +  gene->getTypeBiotype()
        + " " + ((gene->getHavanaTypeId().size() > 0) ? gene->getHavanaTypeId() : "-") + ")";
}

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const Feature* feature) const {
    return fGenomeTransMap->haveQuerySeq(feature->getSeqid());
}

void GeneMapper::debugRecordMapped(const Feature* feature,
                                   const string& desc,
                                   const string& key) const {
    if (gVerbose) {
        cerr << desc << " " << featureDesc(feature);
        if (key != "") {
            cerr << " key: " << key;
        }
        cerr << endl;
    }
}

/* record gene and it's transcripts as being mapped */
void GeneMapper::recordGeneMapped(const Feature* gene) {
    assert(gene->isGene());
    fMappedIdsNames.insert(getBaseId(gene->getTypeId()));
    debugRecordMapped(gene, "recordGeneMapped typeId", getBaseId(gene->getTypeId()));
    
    // N.B.  Don't use gene name for automatic non-coding, as some small
    // non-coding genes has the same name for multiple instances
    // N.B. gene names with `.' are not always a version
    if (not gene->isAutomaticSmallNonCodingGene()) {
        fMappedIdsNames.insert(gene->getTypeName());
        debugRecordMapped(gene, "recordGeneMapped typeName", getBaseId(gene->getTypeName()));
    }
    if (gene->getHavanaTypeId() != "") {
        fMappedIdsNames.insert(getBaseId(gene->getHavanaTypeId()));
        debugRecordMapped(gene, "recordGeneMapped havanaTypeId", getBaseId(gene->getHavanaTypeId()));
    }

    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        recordTranscriptMapped(gene->getChild(i));
    }
}


/* record transcript as being mapped */
void GeneMapper::recordTranscriptMapped(const Feature* transcript) {
    assert(transcript->isTranscript());
    fMappedIdsNames.insert(getBaseId(transcript->getTypeId()));
    debugRecordMapped(transcript, "recordTranscriptMapped typeId", getBaseId(transcript->getTypeId()));
    if (transcript->getHavanaTypeId() != "") {
        fMappedIdsNames.insert(getBaseId(transcript->getHavanaTypeId()));
        debugRecordMapped(transcript, "recordTranscriptMapped havanaTypeId", getBaseId(transcript->getHavanaTypeId()));
    }
}


/* check if gene have been mapped */
bool GeneMapper::checkGeneMapped(const Feature* gene) const {
    assert(gene->isGene());
    if (fMappedIdsNames.find(getBaseId(gene->getTypeId())) != fMappedIdsNames.end()) {
        debugRecordMapped(gene, "checkGeneMapped found typeId", getBaseId(gene->getTypeId()));
        return true;
    }
    if ((not gene->isAutomaticSmallNonCodingGene())
        and fMappedIdsNames.find(getBaseId(gene->getTypeName())) != fMappedIdsNames.end()) {
        debugRecordMapped(gene, "checkGeneMapped found typeName", gene->getTypeName());
        return true;

    }
    if (gene->getHavanaTypeId() != "") {
        if (fMappedIdsNames.find(getBaseId(gene->getHavanaTypeId())) != fMappedIdsNames.end()) {
            debugRecordMapped(gene, "checkGeneMapped found havanaTypeId", gene->getHavanaTypeId());
            return true;
        }
    }
    debugRecordMapped(gene, "checkGeneMapped not found", gene->getTypeId());
    return false;
}

/* check if transcript have already been mapped */
bool GeneMapper::checkTranscriptMapped(const Feature* transcript) const {
    assert(transcript->isTranscript());
    if (fMappedIdsNames.find(getBaseId(transcript->getTypeId())) != fMappedIdsNames.end()) {
        debugRecordMapped(transcript, "checkTranscriptMapped found typeId", getBaseId(transcript->getTypeId()));
        return true;
    }
    return false;
}

/* check if all transcripts in a gene have already been mapped (doesn't check gene) */
bool GeneMapper::checkGeneTranscriptsMapped(const Feature* gene) const {
    assert(gene->isGene());
    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        if (not checkTranscriptMapped(gene->getChild(i))) {
            debugRecordMapped(gene, "checkGeneTranscriptsMapped some transcripts not mapped");
            return false;
        }
    }
    debugRecordMapped(gene, "checkGeneTranscriptsMapped all transcripts mapped");
    return true;
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
        if (transcript->getType() != Feature::TRANSCRIPT) {
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
        cerr << "forceToUnmappedDueToRemapStatus " << featureDesc(mappedGene->src)
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
        cerr << "forceToUnmappedDueToTargetStatus " << featureDesc(mappedGene->src)
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
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " missing target gene due to gene id rename" << endl;
        }
        return false;
    }
    if (not isSrcSeqInMapping(targetGene)) {
        return false;  // sequence not being mapped (moved chroms)
    }

    if (checkForPathologicalGeneRename(mappedGene, targetGene)) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " pathological gene rename" << endl;
        }
        return false;
    }
    
    // allow for pseudogene biotype compatiblity between less and more
    // specific pseudogene biotypes
    if ((targetGene->getTypeBiotype() == mappedGene->src->getTypeBiotype())
        or (targetGene->isPseudogene() == mappedGene->src->isPseudogene())) {
        if (checkGeneTranscriptsMapped(targetGene)) {
            if (gVerbose) {
                cerr << "shouldSubstituteTarget: all transcripts already mapped for target: " << featureDesc(targetGene) << endl;
            }
            return false;
        }
        
#if 0
        // FIXME: This is disables to due to gene ids being reused. OTT gene
        // id reused for DUX4L1, Ensemble gene ids changed too.  This causes
        // both failure to map due to non-overlap and failure to substitute
        // because of thinking it's already mapped.   See idChangeV25Test case.

        
        // must check to make sure it wasn't already mapped due to gene id/name pairing incorrectly changing
        bool alreadyMapped = checkGeneMapped(targetGene);
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: " << !alreadyMapped << ": mapped: " << featureDesc(mappedGene->src)
                 << " substitute target: " << featureDesc(targetGene) << " alreadyMapped: " << alreadyMapped << endl;
        }
        return !alreadyMapped;
#else
        return true;
#endif
    } else {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " gene biotype change: target: " << featureDesc(targetGene)
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
    assert(mappedTranscript->getType() == Feature::TRANSCRIPT);
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
        cerr << "buildGeneFeature: "  << featureDesc(srcGeneTree)
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

/* get target annotation  for a feature, if available */
const Feature* GeneMapper::getTargetAnnotation(const Feature* feature) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    }
    // Try id, havana id, then name.  Only havana id and name are check for
    // gene features.  Transcripts can move between genes keeping the same
    // havana id and transcript names are just a numbering withing the gene
    // and not stable.
    const Feature* targetFeature = fTargetAnnotations->getFeatureById(getPreMappedId(feature->getTypeId()),
                                                                      feature->getSeqid());
    if (feature->isGene()) {
        if ((targetFeature == NULL) and (feature->getHavanaTypeId() != "")) {
            targetFeature = fTargetAnnotations->getFeatureByName(getPreMappedId(feature->getHavanaTypeId()),
                                                                 feature->getSeqid());
        }
        if (targetFeature == NULL) {
            targetFeature = fTargetAnnotations->getFeatureByName(feature->getTypeName(),
                                                                 feature->getSeqid());
        }
    }
    return targetFeature;
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
TargetStatus GeneMapper::getTargetAnnotationStatus(const ResultFeatures* mappedFeature) const {
    if (fTargetAnnotations == NULL) {
        return TARGET_STATUS_NA;
    }
    const Feature* targetFeature = getTargetAnnotation(mappedFeature->src);
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
    const Feature* targetFeature = getTargetAnnotation(mappedFeature->src);
    if (targetFeature == NULL) {
        return emptyString;
    } else {
        return targetFeature->getTypeBiotype();
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
    outputSrcGeneInfo(&mappedGene, mappingInfoFh);

    // must be done after forcing status above
    if (mappedGene.mapped == NULL) {
        outputUnmappedGeneInfo(&mappedGene, mappingInfoFh);
        if (shouldSubstituteTarget(&mappedGene)) {
            substituteTarget(&mappedGene);
            outputTargetGeneInfo(&mappedGene, "targetSubst", mappingInfoFh);
        }
    }
    if (mappedGene.mapped != NULL) {
        featureTreePolish.polishGene(mappedGene.mapped);
        outputMappedGeneInfo(&mappedGene, mappingInfoFh);  // MUST do before saveGene, as it moved to output sets
    }
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
        cerr << "shouldIncludeTargetGene: " << featureDesc(targetGene)
             << " noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(targetGene))
             << " shouldMapGeneType: " << shouldMapGeneType(targetGene)
             << " already mapped: " << checkGeneMapped(targetGene)
             << " gene transcripts all mapped: " << checkGeneTranscriptsMapped(targetGene)
             << endl;
    }
    if (not shouldMapGeneType(targetGene)) {
        // biotypes not excluding from mapped, checkGeneMapped handles
        // case where biotype has changed
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: isSrcSeqInMapping: " << isSrcSeqInMapping(targetGene) << endl;
        }
        return isSrcSeqInMapping(targetGene) and (not checkGeneMapped(targetGene)) and (not checkGeneTranscriptsMapped(targetGene));
    }
    if ((fUseTargetFlags & useTargetForPatchRegions) && inTargetPatchRegion(targetGene)) {
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: in patched region: "
                 << " overlaps mapping: " << checkTargetOverlappingMapped(targetGene, mappedSet) << endl;
        }
        // don't use if there is a mapped with significant overlap
        return (not checkGeneMapped(targetGene)) and (not checkGeneTranscriptsMapped(targetGene))
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
    if (gVerbose) {
        cerr << "copyTargetGene " << featureDesc(targetGene) << endl;
    }
    ResultFeatures mappedGene;
    mappedGene.target = targetGene->cloneTree();
    mappedGene.target->rsetRemapStatus(getNoMapRemapStatus(targetGene));
    mappedGene.target->rsetRemapStatusAttr();
    mappedGene.target->rsetTargetStatusAttr();
    mappedGene.target->rsetSubstitutedMissingTargetAttr(fSubstituteTargetVersion);
    outputTargetGeneInfo(&mappedGene, "targetCopy", mappingInfoFh); // MUST do before copying
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
            fCurrentGeneNum++;
            copyTargetGene(genes[iGene], mappedSet, mappingInfoFh);
        }
    }
}

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfWriter& mappedGxfFh,
                        GxfWriter* unmappedGxfFh,
                        ostream& mappingInfoFh,
                        ostream* transcriptPslFh) {
    AnnotationSet mappedSet(&fGenomeTransMap->fTargetSizes);
    AnnotationSet unmappedSet(&fGenomeTransMap->fQuerySizes);
    FeatureTreePolish featureTreePolish(fPreviousMappedAnotations);
    
    const FeatureVector& srcGenes = fSrcAnnotations->getGenes();
    outputInfoHeader(mappingInfoFh);
    for (int i = 0; i < srcGenes.size(); i++) {
        if (gVerbose) {
            cerr << endl << "mapGxf: " << featureDesc(srcGenes[i])
                 << " shouldMapGeneType: " << shouldMapGeneType(srcGenes[i])
                 << " noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(srcGenes[i]))
                 << " " << srcGenes[i]->getTypeId() << " " << srcGenes[i]->getSource()
                 << endl;
        }
        if (shouldMapGeneType(srcGenes[i])) {
            fCurrentGeneNum++;
            mapGene(srcGenes[i], mappedSet, unmappedSet, featureTreePolish, mappingInfoFh, transcriptPslFh);
        }
    }
    if ((fUseTargetFlags != 0) and (fTargetAnnotations != NULL)) {
        copyTargetGenes(mappedSet, mappingInfoFh);
    }
    mappedSet.sort();
    mappedSet.write(mappedGxfFh);
    if (unmappedGxfFh != NULL) {
        unmappedSet.sort();
        unmappedSet.write(*unmappedGxfFh);
    }
}

