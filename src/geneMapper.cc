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
const float geneExpansionThreshold = 0.20;  

/*  mapinfo TSV headers, terminated by NULL */
static const char* mappingInfoHeaders[] = {
    "id", "name", "type", "biotype", "source",
    "srcChrom", "srcStart", "srcEnd", "srcStrand",
    "mappedChrom", "mappedStart", "mappedEnd", "mappedStrand",
    "mappingStatus", "numMappings",
    "targetStatus", "targetBiotype", NULL
};


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
void GeneMapper::forceToUnmapped(FeatureNode* featureNode,
                                 RemapStatus remapStatus) const {
    FeatureMapper::forceToUnmapped(featureNode);
    featureNode->setRemapStatus(remapStatus);
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        forceToUnmapped(featureNode->fChildren[i], remapStatus);
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
    return (start == 0) ? 0 : (end-start)+1;
}

/* check if gene transcripts cause the gene size to expand beyond a
 * threshold*/
bool GeneMapper::hasExcessiveSizeChange(FeatureNode* geneTree) const {
    int srcGeneLength = geneTree->fFeature->size();
    int mappedGeneLength = calcMappedGeneLength(geneTree);
    return (abs(float(mappedGeneLength-srcGeneLength)) / float(srcGeneLength)) > geneExpansionThreshold;
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
    FeatureMapper::mapBounding(geneTree, seqid, start, end, strand);
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
}

/* output GFF3 mapped ##sequence-region if not already written */
void GeneMapper::outputMappedSeqRegionIfNeed(FeatureNode* geneTree,
                                             ostream& mappedGxfFh) {
    if (geneTree->fFeature->getFormat() == GFF3_FORMAT) {
        const string& mappedSeqId = geneTree->fMappedFeatures[0]->fSeqid;
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
        return fTargetAnnotations->get(featureNode->fFeature->getTypeId(),
                                       featureNode->fFeature->fSeqid);
    }
}

/* If target gene annotations are available, get status of mapping
 * relative to older version of gene. */
const string& GeneMapper::getTargetAnnotationStatus(FeatureNode* featureNode,
                                                    const GxfFeature* targetFeature) const {
    static const string STATUS_NA = "na";
    static const string STATUS_NEW = "new";
    static const string STATUS_LOST = "lost";
    static const string STATUS_OVERLAP = "overlap";
    static const string STATUS_NONOVERLAP = "nonOverlap";
    if (fTargetAnnotations == NULL) {
        return STATUS_NA;
    }
    if (targetFeature == NULL) {
        return STATUS_NEW;
    }
    if (featureNode->fMappedFeatures.size() == 0) {
        return STATUS_LOST;
    }
    assert(featureNode->fMappedFeatures.size() == 1);  // genes/transcripts don't split
    if (featureNode->fMappedFeatures[0]->overlaps(targetFeature)) {
        return STATUS_OVERLAP;
    } else {
        return STATUS_NONOVERLAP;
    }
}

/* If target gene annotations are available, get biotype of target feature */
const string& GeneMapper::getTargetAnnotationBiotype(FeatureNode* featureNode,
                                                     const GxfFeature* targetFeature) const {
    if (targetFeature == NULL) {
        return emptyString;
    } else {
        return targetFeature->getTypeBiotype();
    }
}

/* output info on one bounding feature */
void GeneMapper::outputFeatureInfo(FeatureNode* featureNode,
                                   ostream& mappingInfoFh) const {
    assert(featureNode->fMappedFeatures.size() <= 1);  // only handles bounding features
    // get target if available
    const GxfFeature* targetFeature = getTargetAnnotation(featureNode);
    
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
                  << getTargetAnnotationStatus(featureNode, targetFeature) << "\t"
                  << getTargetAnnotationBiotype(featureNode, targetFeature) << endl;
}

/*
 * Output information about gene mapping
 */
void GeneMapper::outputInfo(FeatureNode* geneNode,
                            ostream& mappingInfoFh) const {
    outputFeatureInfo(geneNode, mappingInfoFh);
    // transcripts
    for (int i = 0; i < geneNode->fChildren.size(); i++) {
        outputFeatureInfo(geneNode->fChildren[i], mappingInfoFh);
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
    buildGeneFeature(geneTree);
    // handle gene level conflicted */
    if (hasMixedMappedSeqStrand(geneTree)) {
        forceToUnmapped(geneTree, REMAP_STATUS_GENE_CONFLICT);
    } else if (hasExcessiveSizeChange(geneTree)) {
        forceToUnmapped(geneTree, REMAP_STATUS_GENE_SIZE_CHANGE);
    } else {
        geneTree->setRemapStatusFromChildren();
    }
    geneTree->setRemapStatusAttr();
    geneTree->setNumMappingsAttr();
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

