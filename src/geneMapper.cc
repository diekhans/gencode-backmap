#include "geneMapper.hh"
#include "featureMapper.hh"
#include "featureTree.hh"
#include "bedMap.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "transcriptMapper.hh"
#include "annotationSet.hh"
#include "featureTreePolish.hh"
#include "globals.hh"
#include "gxf.hh"


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
                            const FeatureNode* feature,
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
 * Output information about source gene that was mapped or failed mapping.
 */
void GeneMapper::outputSrcGeneInfo(const ResultFeatureTrees* mappedGene,
                                   ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    const FeatureNode* srcGene = mappedGene->src;
    outputInfo("mapSrc", "gene", srcGene, srcGene->getRemapStatus(), 0, srcGene->getTargetStatus(), mappingInfoFh);
    for (int i = 0; i < srcGene->getNumChildren(); i++) {
        const FeatureNode* srcTrans = srcGene->getChild(i);
        outputInfo("mapSrc", "trans", srcTrans, srcTrans->getRemapStatus(), 0, srcTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene that was mapped or partially mapped.
 */
void GeneMapper::outputMappedGeneInfo(const ResultFeatureTrees* mappedGene,
                                      ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    // not all transcripts maybe not be mapped
    const FeatureNode* mapGene = mappedGene->mapped;
    outputInfo("map", "gene", mapGene, mapGene->getRemapStatus(), mapGene->getNumMappings(), mapGene->getTargetStatus(), mappingInfoFh);

    // transcripts
    for (int i = 0; i < mapGene->getNumChildren(); i++) {
        const FeatureNode* mapTrans = mapGene->getChild(i);
        outputInfo("map", "trans", mapTrans, mapTrans->getRemapStatus(), mapGene->getNumMappings(), mapTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene that was unmapped or partially unmapped.
 */
void GeneMapper::outputUnmappedGeneInfo(const ResultFeatureTrees* mappedGene,
                                        ostream& mappingInfoFh) const {
    assert(mappedGene->src != NULL);
    // not all transcripts maybe not be mapped
    const FeatureNode* unmapGene = mappedGene->unmapped;
    outputInfo("unmap", "gene", unmapGene, unmapGene->getRemapStatus(), unmapGene->getNumMappings(), unmapGene->getTargetStatus(), mappingInfoFh);

    // transcripts
    for (int i = 0; i < unmapGene->getNumChildren(); i++) {
        const FeatureNode* unmapTrans = unmapGene->getChild(i);
        outputInfo("unmap", "trans", unmapTrans, unmapTrans->getRemapStatus(), unmapTrans->getNumMappings(), unmapTrans->getTargetStatus(), mappingInfoFh);
    }
}

/*
 * Output information about gene where target is used by default or substituted
 */
void GeneMapper::outputTargetGeneInfo(const ResultFeatureTrees* mappedGene,
                                      const string& targetAction,
                                      ostream& mappingInfoFh) const {
    const FeatureNode* targetGene = mappedGene->target;
    outputInfo(targetAction, "gene", targetGene, targetGene->getRemapStatus(), 0, targetGene->getTargetStatus(), mappingInfoFh);
    for (int i = 0; i < targetGene->getNumChildren(); i++) {
        const FeatureNode* targetTrans = targetGene->getChild(i);
        outputInfo(targetAction, "trans", targetTrans, targetTrans->getRemapStatus(), 0, targetTrans->getTargetStatus(), mappingInfoFh);
    }
}

/* get string describing gene or other feature */
string GeneMapper::featureDesc(const FeatureNode* feature) const {
    return "(" + feature->getTypeId() + " " + feature->fFeature->getSeqid() +
        " " + feature->getTypeName() + " " +  feature->getTypeBiotype() +
        " " + ((not feature->getHavanaTypeId().empty()) ? feature->getHavanaTypeId() : "-") + ")";
}

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const FeatureNode* feature) const {
    return fGenomeTransMap->haveQuerySeq(feature->getSeqid());
}

void GeneMapper::debugRecordMapped(const FeatureNode* feature,
                                   const string& desc,
                                   const string& key) const {
    if (gVerbose) {
        cerr << desc << " " << featureDesc(feature);
        if (key != "") {
            cerr << " key: (" << key << ", " << feature->getSeqid() << ")";
        }
        cerr << endl;
    }
}

/* record gene and it's transcripts as being mapped */
void GeneMapper::recordGeneMapped(const FeatureNode* gene) {
    assert(gene->isGene());
    fMappedIdsNames.addBaseId(gene->getTypeId(), gene->fFeature->getSeqid());
    debugRecordMapped(gene, "recordGeneMapped typeId", getBaseId(gene->getTypeId()));
    
    if (useGeneNameForMappingKey(gene)) {
        fMappedIdsNames.addName(gene->getTypeName(), gene->fFeature->getSeqid());
        debugRecordMapped(gene, "recordGeneMapped typeName", gene->getTypeName());
    }
    if (gene->getHavanaTypeId() != "") {
        fMappedIdsNames.addBaseId(gene->getHavanaTypeId(), gene->fFeature->getSeqid());
        debugRecordMapped(gene, "recordGeneMapped havanaTypeId", getBaseId(gene->getHavanaTypeId()));
    }
}


/* remove record of gene and it's transcripts as being mapped */
void GeneMapper::forgetGeneMapped(const FeatureNode* gene) {
    assert(gene->isGene());
    fMappedIdsNames.removeBaseId(gene->getTypeId(), gene->fFeature->getSeqid());
    debugRecordMapped(gene, "forgetGeneMapped typeId", getBaseId(gene->getTypeId()));
    
    if (useGeneNameForMappingKey(gene)) {
        fMappedIdsNames.removeName(gene->getTypeName(), gene->fFeature->getSeqid());
        debugRecordMapped(gene, "forgetGeneMapped typeName", gene->getTypeName());
    }
    if (gene->getHavanaTypeId() != "") {
        fMappedIdsNames.removeBaseId(gene->getHavanaTypeId(), gene->fFeature->getSeqid());
        debugRecordMapped(gene, "forgetGeneMapped havanaTypeId", getBaseId(gene->getHavanaTypeId()));
    }
}


/* record transcript as being mapped */
void GeneMapper::recordTranscriptMapped(const FeatureNode* transcript) {
    assert(transcript->isTranscript());
    if (fMappedIdsNames.haveBaseId(transcript->getTypeId(), transcript->fFeature->getSeqid())) {
        //FIXME: throw logic_error("BUG: transcript base id has already been mapped: " + transcript->getTypeId());
        cerr  << "WARNING: transcript base id has already been mapped: " + transcript->getTypeId() << endl;
    }
    fMappedIdsNames.addBaseId(transcript->getTypeId(), transcript->fFeature->getSeqid());
    debugRecordMapped(transcript, "recordTranscriptMapped typeId", getBaseId(transcript->getTypeId()));
    if (transcript->getHavanaTypeId() != "") {
        fMappedIdsNames.addBaseId(transcript->getHavanaTypeId(), transcript->fFeature->getSeqid());
        debugRecordMapped(transcript, "recordTranscriptMapped havanaTypeId", getBaseId(transcript->getHavanaTypeId()));
    }
}

/* removed record of transcript as being mapped */
void GeneMapper::forgetTranscriptMapped(const FeatureNode* transcript) {
    assert(transcript->isTranscript());
    fMappedIdsNames.removeBaseId(transcript->getTypeId(), transcript->fFeature->getSeqid());
    debugRecordMapped(transcript, "forgetTranscriptMapped typeId", getBaseId(transcript->getTypeId()));
    if (transcript->getHavanaTypeId() != "") {
        fMappedIdsNames.removeBaseId(transcript->getHavanaTypeId(), transcript->fFeature->getSeqid());
        debugRecordMapped(transcript, "forgetTranscriptMapped havanaTypeId", getBaseId(transcript->getHavanaTypeId()));
    }
}


/* check if gene have been mapped */
bool GeneMapper::checkGeneMapped(const FeatureNode* gene) const {
    assert(gene->isGene());
    if (fMappedIdsNames.haveBaseId(gene->getTypeId(), gene->fFeature->getSeqid())) {
        debugRecordMapped(gene, "checkGeneMapped found typeId", getBaseId(gene->getTypeId()));
        return true;
    }
    if (useGeneNameForMappingKey(gene)
        and fMappedIdsNames.haveName(gene->getTypeName(), gene->fFeature->getSeqid())) {
        debugRecordMapped(gene, "checkGeneMapped found typeName", gene->getTypeName());
        return true;

    }
    if (gene->getHavanaTypeId() != "") {
        if (fMappedIdsNames.haveBaseId(gene->getHavanaTypeId(), gene->fFeature->getSeqid())) {
            debugRecordMapped(gene, "checkGeneMapped found havanaTypeId", gene->getHavanaTypeId());
            return true;
        }
    }
    debugRecordMapped(gene, "checkGeneMapped not found", gene->getTypeId());
    return false;
}

/* check if transcript have already been mapped */
bool GeneMapper::checkTranscriptMapped(const FeatureNode* transcript) const {
    assert(transcript->isTranscript());
    bool found = fMappedIdsNames.haveBaseId(transcript->getTypeId(), transcript->fFeature->getSeqid());
    debugRecordMapped(transcript, std::string("checkTranscriptMapped ") + (found ? "TRUE" : "FALSE")  + " typeId", getBaseId(transcript->getTypeId()));
    return found;
}

/* check if all transcripts in a gene have already been mapped (doesn't check gene) */
bool GeneMapper::checkAllGeneTranscriptsMapped(const FeatureNode* gene) const {
    assert(gene->isGene());
    for (size_t i = 0; i < gene->getNumChildren(); i++) {
        if (not checkTranscriptMapped(gene->getChild(i))) {
            debugRecordMapped(gene, "checkAllGeneTranscriptsMapped some transcripts not mapped");
            return false;
        }
    }
    debugRecordMapped(gene, "checkAllGeneTranscriptsMapped all transcripts mapped");
    return true;
}

/* check if any transcripts in a gene have already been mapped (doesn't check gene) */
bool GeneMapper::checkAnyGeneTranscriptsMapped(const FeatureNode* gene) const {
    assert(gene->isGene());
    for (size_t i = 0; i < gene->getNumChildren(); i++) {
        if (checkTranscriptMapped(gene->getChild(i))) {
            debugRecordMapped(gene, "checkAnyGeneTranscriptsMapped some transcripts mapped");
            return true;
        }
    }
    debugRecordMapped(gene, "checkAnyGeneTranscriptsMapped no transcripts mapped");
    return false;
}

/* process one transcript */
ResultFeatureTrees GeneMapper::processTranscript(const FeatureNode* transcript,
                                                 ostream* transcriptPslFh) {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcript, fTargetAnnotations,
                                      isSrcSeqInMapping(transcript), transcriptPslFh);
    ResultFeatureTrees mappedTranscript = transcriptMapper.mapTranscriptFeatures(transcript);
    TargetStatus targetStatus = getTargetAnnotationStatus(&mappedTranscript);
    mappedTranscript.setTargetStatus(targetStatus);
    if (mappedTranscript.mapped != NULL) {
        recordTranscriptMapped(transcript);
    }
    return mappedTranscript;
}

/* process all transcripts of gene. */
ResultFeatureTreesVector GeneMapper::processTranscripts(const FeatureNode* gene,
                                                        ostream* transcriptPslFh) {
    ResultFeatureTreesVector mappedTranscripts;
    for (size_t i = 0; i < gene->getNumChildren(); i++) {
        const FeatureNode* transcript = gene->getChild(i);
        if (transcript->getType() != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcript->toString());
        }
        mappedTranscripts.push_back(processTranscript(transcript, transcriptPslFh));
    }
    return mappedTranscripts;
}

/* find a matching gene or transcript given by id */
FeatureNode* GeneMapper::findMatchingBoundingFeature(const FeatureNodeVector& features,
                                                     const FeatureNode* feature) const {
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
void GeneMapper::copyMappingMetadata(const FeatureNode* origFeature,
                                     FeatureNode* newFeature) const {
    newFeature->setRemapStatus(origFeature->getRemapStatus());
    newFeature->setTargetStatus(origFeature->getTargetStatus());
    newFeature->setNumMappings(max(origFeature->getNumMappings(), newFeature->getNumMappings()));
}

/* copy gene and transcript attributes to a fresh unmapped clone,
 * only at the transcript and gene levels. */
void GeneMapper::copyGeneMetadata(const FeatureNode* origGene,
                                  FeatureNode* newGene) const {
    copyMappingMetadata(origGene, newGene);

    // copy transcript metadata
    for (size_t i = 0; i < newGene->getNumChildren(); i++) {
        FeatureNode* newTranscript = newGene->getChild(i);
        const FeatureNode* origTranscript = findMatchingBoundingFeature(origGene->getChildren(), newTranscript);
        if (origTranscript != NULL) {
            copyMappingMetadata(origTranscript, newTranscript);
        }
    }
}

/* force to unmapped */
void GeneMapper::forceToUnmapped(ResultFeatureTrees* mappedGene) {
    // create copy and update attributes at gene/transcript level before freeing.
    FeatureNode* srcCopy = mappedGene->src->cloneTree();
    if (mappedGene->mapped != NULL) {
        copyGeneMetadata(mappedGene->mapped, srcCopy);
        for (size_t i = 0; i < mappedGene->mapped->getNumChildren(); i++) {
            forgetTranscriptMapped(mappedGene->mapped->getChild(i));
        }
        forgetGeneMapped(mappedGene->mapped);
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
                                                 RemapStatus remapStatus) {
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
void GeneMapper::forceToUnmappedDueToTargetStatus(ResultFeatureTrees* mappedGene,
                                                  TargetStatus targetStatus) {
    if (gVerbose) {
        cerr << "forceToUnmappedDueToTargetStatus " << featureDesc(mappedGene->src)
             << " " << targetStatusToStr(targetStatus) << endl;
    }
    forceToUnmapped(mappedGene);
    mappedGene->setTargetStatus(targetStatus);
    mappedGene->rsetTargetStatusAttr();
}

/* check if gene contains transcripts with different mapped seqid/strand  */
bool GeneMapper::hasMixedMappedSeqStrand(const ResultFeatureTrees* mappedGene) const {
    int numTranscripts = (mappedGene->mapped != NULL) ? mappedGene->mapped->getNumChildren() : 0;
    string seqid, strand;
    for (int i = 0; i < numTranscripts; i++) {
        const FeatureNode* transcript = mappedGene->mapped->getChild(i);
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
bool GeneMapper::hasTargetStatusNonOverlap(const ResultFeatureTrees* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return false;
    }
    if (mappedGene->mapped->getTargetStatus() == TARGET_STATUS_NONOVERLAP) {
        // handles weird case of ENSG00000239810 were gene doesn't overlap
        // however the only transcript is `new'.
        return true;
    }
    for (int i = 0; i < mappedGene->mapped->getNumChildren(); i++) {
        if (mappedGene->mapped->getChild(i)->getTargetStatus() == TARGET_STATUS_NONOVERLAP) {
            return true;
        }
    }
    return false;
}

/* check if gene transcripts cause the gene size to expand beyond a
 * threshold*/
bool GeneMapper::hasExcessiveSizeChange(const ResultFeatureTrees* mappedGene) const {
    if (mappedGene->mapped == NULL) {
        return false;  // not mapped
    } else {
        int srcGeneLength = mappedGene->src->length();
        int mappedGeneLength = mappedGene->mapped->length();
        return (abs(float(mappedGeneLength-srcGeneLength)) / float(srcGeneLength)) > geneExpansionThreshold;
    }
}

/* set the number of mappings of the gene  */
void GeneMapper::setNumGeneMappings(FeatureNode* mappedGeneTree) const {
    for (int i = 0; i < mappedGeneTree->getNumChildren(); i++) {
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
bool GeneMapper::checkForPathologicalGeneRename(const ResultFeatureTrees* mappedGene,
                                                const FeatureNode* targetGene) const {
    return (getBaseId(mappedGene->src->getTypeId()) != getBaseId(targetGene->getTypeId()))
        and (fSrcAnnotations->getFeatureByIdChrom(targetGene->getTypeId(), targetGene->fFeature->getSeqid()) != NULL);
}


/* should we substitute target version of gene even if mapped?  */
bool GeneMapper::shouldSubstituteTarget(const ResultFeatureTrees* mappedGene) const {
    // mis-mapped gene or the right biotype to a mapped sequence
    const FeatureNode* targetGene = getTargetAnnotation(mappedGene->src);
    if (targetGene == NULL) {
        // this should not have happened, but it does because of a transcript id
        // ENST00000426406 incorrectly being moved to a different gene
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " missing target gene probably due to gene id rename" << endl;
        }
        return false;
    }
    if (not isSrcSeqInMapping(targetGene)) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " target gene chrom not in mapping alignments" << endl;
        }
        return false;  // sequence not being mapped (moved chroms)
    }

    if (checkForPathologicalGeneRename(mappedGene, targetGene)) {
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                 << " pathological gene rename" << endl;
        }
        return false;
    }

    // FIXME: tmp hack for V31 to avoid substitution of target when gene has merged
    // this is being fixed properly , but needed quickly
    if (targetGene->getTypeId() == "ENSG00000226155.1") {
#if 0
        cerr << "NOTE: V31 hack, don't include target gene "  << targetGene->getTypeId() << endl;
        return false;
#else
        cerr << "NOTE: V31 hack DISABLED "  << targetGene->getTypeId() << endl;
#endif
    }
    
    // allow for pseudogene biotype compatiblity between less and more
    // specific pseudogene biotypes
    if ((targetGene->getTypeBiotype() == mappedGene->src->getTypeBiotype())
        or (targetGene->isPseudogene() == mappedGene->src->isPseudogene())) {
        if (checkAllGeneTranscriptsMapped(targetGene)) {
            if (gVerbose) {
                cerr << "shouldSubstituteTarget: false: " << featureDesc(mappedGene->src)
                     << ", some transcripts already mapped for target" << featureDesc(targetGene) << endl;
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
        if (gVerbose) {
            cerr << "shouldSubstituteTarget: true: " << featureDesc(mappedGene->src)
                 << " with: " << featureDesc(targetGene) << endl;
        }
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
void GeneMapper::substituteTarget(ResultFeatureTrees* mappedGene) {
    if (gVerbose) {
        cerr << "substituteTarget: " << featureDesc(mappedGene->src)
             << " with: " << featureDesc(getTargetAnnotation(mappedGene->src))
             << endl;
    }
    assert(mappedGene->target == NULL);
    mappedGene->target = getTargetAnnotation(mappedGene->src)->cloneTree();
    for (int i = 0; i < mappedGene->target->getNumChildren(); i++) {
        recordTranscriptMapped(mappedGene->target->getChild(i));
    }

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
        assert(mappedGene.mapped == NULL);
        mappedGene.mapped = buildMappedGeneFeature(srcGeneTree, mappedTranscripts);
        
    }
    if (mappedTranscripts.haveUnmapped()) {
        assert(mappedGene.unmapped == NULL);
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
void GeneMapper::saveMapped(ResultFeatureTrees& mappedGene,
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
void GeneMapper::saveUnmapped(ResultFeatureTrees& mappedGene,
                              AnnotationSet& unmappedSet) {
    if (mappedGene.unmapped) {
        unmappedSet.addGene(mappedGene.unmapped);
        mappedGene.unmapped = NULL;
    }
}

/* get target annotation  for a feature, if available */
const FeatureNode* GeneMapper::getTargetAnnotation(const FeatureNode* feature) const {
    if (fTargetAnnotations == NULL) {
        return NULL;
    }
    // Try id, havana id, then name.  Only havana id and name are check for
    // gene features.  Transcripts can move between genes keeping the same
    // havana id and transcript names are just a numbering withing the gene
    // and not stable.
    const FeatureNode* targetFeature = fTargetAnnotations->getFeatureByIdChrom(getPreMappedId(feature->getTypeId()),
                                                                               feature->fFeature->getSeqid());
    if (feature->isGene()) {
        if ((targetFeature == NULL) and (feature->getHavanaTypeId() != "")) {
            targetFeature = fTargetAnnotations->getFeatureByNameChrom(getPreMappedId(feature->getHavanaTypeId()),
                                                                      feature->fFeature->getSeqid());
        }
        if ((targetFeature == NULL) and useGeneNameForMappingKey(feature)) {
            targetFeature = fTargetAnnotations->getFeatureByNameChrom(feature->getTypeName(),
                                                                      feature->fFeature->getSeqid());
        }
    }
    return targetFeature;
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
TargetStatus GeneMapper::getTargetAnnotationStatus(const ResultFeatureTrees* mappedFeature) const {
    if (fTargetAnnotations == NULL) {
        return TARGET_STATUS_NA;
    }
    const FeatureNode* targetFeature = getTargetAnnotation(mappedFeature->src);
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
const string& GeneMapper::getTargetAnnotationBiotype(const ResultFeatureTrees* mappedFeature) const {
    const FeatureNode* targetFeature = getTargetAnnotation(mappedFeature->src);
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
void GeneMapper::processGeneLevelMapping(ResultFeatureTrees* mappedGene) {
    if (hasMixedMappedSeqStrand(mappedGene)) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_GENE_CONFLICT);
    } else if (hasExcessiveSizeChange(mappedGene)) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_GENE_SIZE_CHANGE);
    } else if (mappedGene->getRemapStatus() == REMAP_STATUS_NON_PRIMARY) {
        forceToUnmappedDueToRemapStatus(mappedGene, REMAP_STATUS_NON_PRIMARY);
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
                         FeatureTreePolish& featureTreePolish,
                         ostream& mappingInfoFh,
                         ostream* transcriptPslFh) {
    if (gVerbose) {
        cerr << "mapGene: "  << featureDesc(srcGeneTree) << endl;
    }
    
    ResultFeatureTreesVector mappedTranscripts = processTranscripts(srcGeneTree, transcriptPslFh);
    ResultFeatureTrees mappedGene = buildGeneFeature(srcGeneTree, mappedTranscripts);
    setGeneLevelMappingAttributes(&mappedGene);
    processGeneLevelMapping(&mappedGene);
    outputSrcGeneInfo(&mappedGene, mappingInfoFh);

    // must be done after forcing status above
    if (mappedGene.mapped == NULL) {
        outputUnmappedGeneInfo(&mappedGene, mappingInfoFh);
        if ((not fSubstituteTargetVersion.empty()) and shouldSubstituteTarget(&mappedGene)) {
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

/*
 * map and output one gene's annotations
 */
void GeneMapper::maybeMapGene(const FeatureNode* srcGeneTree,
                              AnnotationSet& mappedSet,
                              AnnotationSet& unmappedSet,
                              FeatureTreePolish& featureTreePolish,
                              ostream& mappingInfoFh,
                              ostream* transcriptPslFh) {
    if (gVerbose) {
        cerr << endl << "maybeMapGene: " << featureDesc(srcGeneTree)
             << " shouldMapGeneType: " << shouldMapGeneType(srcGeneTree)
             << " noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(srcGeneTree))
             << " " << srcGeneTree->getTypeId() << " " << srcGeneTree->getSource()
             << endl;
    }
    if (shouldMapGeneType(srcGeneTree)) {
        if (checkAnyGeneTranscriptsMapped(srcGeneTree)) {
            // FIXME: should output special status
            cerr << "WARNING: " << srcGeneTree->getTypeId() << " has transcript mapped in other genes" << endl;
            outputInfo("mapSrc", "gene", srcGeneTree, REMAP_STATUS_ERROR, 0, TARGET_STATUS_ERROR, mappingInfoFh);
            for (int i = 0; i < srcGeneTree->getNumChildren(); i++) {
                outputInfo("mapSrc", "trans", srcGeneTree->getChild(i), REMAP_STATUS_ERROR, 0, TARGET_STATUS_ERROR, mappingInfoFh);
            }
        } else {
            fCurrentGeneNum++;
            mapGene(srcGeneTree, mappedSet, unmappedSet, featureTreePolish, mappingInfoFh, transcriptPslFh);
        }
    }
    if (gVerbose) {
        cerr << "----" << endl;
    }
}

/* determine if this is a gene type that should not be mapped, returning
 * the remap status */
RemapStatus GeneMapper::getNoMapRemapStatus(const FeatureNode* gene) const {
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
bool GeneMapper::shouldMapGeneType(const FeatureNode* gene) const {
    return getNoMapRemapStatus(gene) == REMAP_STATUS_NONE;
}

/* is target gene in patch region? */
bool GeneMapper::inTargetPatchRegion(const FeatureNode* targetGene) {
    return fTargetPatchMap->anyOverlap(targetGene->getSeqid(),
                                       targetGene->getStart(),
                                       targetGene->getEnd());
}

/* check to see if the target overlaps a mapped gene with sufficient similarity to
 * be considered the same annotation..  */
bool GeneMapper::checkTargetOverlappingMapped(const FeatureNode* targetGene,
                                              AnnotationSet& mappedSet) {
    static const float minSimilarity = 0.5;
    FeatureNodeVector overlapping = mappedSet.findOverlappingGenes(targetGene, minSimilarity,
                                                                   fOnlyManualForTargetSubstituteOverlap);
    return not overlapping.empty();
}

/* check for special cases where target should not be included */
bool GeneMapper::checkForIncludeTargetSpecialCases(const FeatureNode* targetGene) const {
    // ENSG00000168939 (SPRY3) was on PAR X/Y in GRCh37, sometime it GRCh38 it
    // was removed from chrY since it had one transcript with an exon outside
    // of the PAR.  So with the the same id for PAR approach,
    // 
    if ((getBaseId(targetGene->getTypeId()) == "ENSG00000168939") and (targetGene->getSeqid() == "chrY")) {
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: false: " << featureDesc(targetGene)
                 << "PAR-Y pathological case" << endl;
        }
        return true;
    }
    return false;
}


/*
 * Check if a target gene should be copied.
 */
bool GeneMapper::shouldIncludeTargetGene(const FeatureNode* targetGene,
                                         AnnotationSet& mappedSet)  {
    if (gVerbose) {
        cerr << "shouldIncludeTargetGene: " << featureDesc(targetGene)
             << ": noMapRemapStatus: " << remapStatusToStr(getNoMapRemapStatus(targetGene))
             << ", shouldMapGeneType: " << shouldMapGeneType(targetGene)
             << ", alreadyMapped: " << checkGeneMapped(targetGene)
             << ", anyTransMapped: " << checkAnyGeneTranscriptsMapped(targetGene) << endl;
    }
    if (checkForIncludeTargetSpecialCases(targetGene)) {
        return false;
    }
    
    if (not shouldMapGeneType(targetGene)) {
        // biotypes not excluding from mapped, checkGeneMapped handles
        // case where biotype has changed
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: isSrcSeqInMapping: " << isSrcSeqInMapping(targetGene)
                 << ", not-checkGeneMapped: " << (not checkGeneMapped(targetGene))
                 << ", not-checkAnyGeneTranscriptsMapped: " << (not checkAnyGeneTranscriptsMapped(targetGene))
                 << endl;
        }
        return isSrcSeqInMapping(targetGene) and (not checkGeneMapped(targetGene)) and (not checkAnyGeneTranscriptsMapped(targetGene));
    }
    if ((fUseTargetFlags & useTargetForPatchRegions) && inTargetPatchRegion(targetGene)) {
        if (gVerbose) {
            cerr << "    shouldIncludeTargetGene: in patched region: "
                 << "not-checkGeneMapped: " << (not checkGeneMapped(targetGene))
                 << ", not checkAnyGeneTranscriptsMapped: " << (not checkAnyGeneTranscriptsMapped(targetGene))
                 << ", checkTargetOverlappingMapped: " << (not checkTargetOverlappingMapped(targetGene, mappedSet)) << endl;

        }
        // don't use if there is a mapped with significant overlap or any of it's transcripts
        // have been mapped
        return (not checkGeneMapped(targetGene)) and (not checkAnyGeneTranscriptsMapped(targetGene))
            and (not checkTargetOverlappingMapped(targetGene, mappedSet));
    }
    if (gVerbose) {
        cerr << "    shouldIncludeTargetGene: FALSE" << endl;
    }
    return false;
}

/*
 * copy a target gene annotation that was skipped for mapping
 */
void GeneMapper::copyTargetGene(const FeatureNode* targetGene,
                                AnnotationSet& mappedSet,
                                ostream& mappingInfoFh) {
    if (gVerbose) {
        cerr << "copyTargetGene " << featureDesc(targetGene) << endl;
    }
    ResultFeatureTrees mappedGene;
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
 * copy a target gene annotation if it should be including.
 */
void GeneMapper::maybeCopyTargetGene(const FeatureNode* targetGene,
                                     AnnotationSet& mappedSet,
                                     ostream& mappingInfoFh) {
    if (gVerbose) {
        cerr << endl << "maybeCopyTargetGene: " << featureDesc(targetGene) << endl;
    }
    if (shouldIncludeTargetGene(targetGene, mappedSet)) {
        fCurrentGeneNum++;
        copyTargetGene(targetGene, mappedSet, mappingInfoFh);
    }
    if (gVerbose) {
        cerr << "----" << endl;
    }
}

/*
 * copy target annotations that are skipped for mapping
 */
void GeneMapper::copyTargetGenes(AnnotationSet& mappedSet,
                                 ostream& mappingInfoFh) {
    const FeatureNodeVector& genes = fTargetAnnotations->getGenes();
    for (int iGene = 0; iGene < genes.size(); iGene++) {
        maybeCopyTargetGene(genes[iGene], mappedSet, mappingInfoFh);
    }
}

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfWriter& mappedGxfFh,
                        ostream& mappingInfoFh,
                        ostream* transcriptPslFh) {
    AnnotationSet mappedSet(&fGenomeTransMap->fTargetSizes);
    AnnotationSet unmappedSet(&fGenomeTransMap->fQuerySizes);
    FeatureTreePolish featureTreePolish(fPreviousMappedAnotations);
    
    const FeatureNodeVector& srcGenes = fSrcAnnotations->getGenes();
    outputInfoHeader(mappingInfoFh);
    for (int i = 0; i < srcGenes.size(); i++) {
        maybeMapGene(srcGenes[i], mappedSet, unmappedSet, featureTreePolish, mappingInfoFh, transcriptPslFh);
    }
    if ((fUseTargetFlags != 0) and (fTargetAnnotations != NULL)) {
        copyTargetGenes(mappedSet, mappingInfoFh);
    }
    mappedSet.sortGencode();
    mappedSet.write(mappedGxfFh);
}

