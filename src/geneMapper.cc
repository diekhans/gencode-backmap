#include "geneMapper.hh"
#include "featureMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "pslOps.hh"
#include "pslMapping.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include "frame.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>

#define debug 0

/* fraction of gene expansion that causes a rejection */
const float geneExpansionThreshold = 0.20;  

/*  TSV headers, terminated by NULL */
static const char* mappingInfoHeaders[] = {
    "id", "type",
    "srcChrom", "srcStart", "srcEnd", "srcStrand",
    "mappedChrom", "mappedStart", "mappedEnd", "mappedStrand",
    "mappingStatus", "numMappings", NULL
};


/**
 * Class to map a single transcript and subfeatures
 * All mapping is done by a two-level alignment.
 * with mapping alignments.
 *  - to map cdsA in  annotation of genomeA to genomeB:
 *    - an alignment of [genomeA->exonsA]
 *    - an alignment of [exonsA->genomeB]
 *    - do a two level transmap:
 *      cdsA->genomeA =>  [genomeA->exonsA] => cdsA->exonA => [exonsA->genomeB] => cdsA->genomeB
 */
class TranscriptMapper {
    private:
    const bool fSrcSeqInMapping;                 // do we have source sequence in genomic mapps
    const PslMapping* fExonsMapping;            // exons as psl and genome mapping of exons.
    TransMapVector fVaiExonsTransMaps;          // transmap objects that are combined
    const FeatureTransMap* fViaExonsFeatureTransMap;   // two-level transmap, NULL if can't map (owned)
    
    /* get exon features */
    static GxfFeatureVector getExons(const FeatureNode* transcriptTree) {
        GxfFeatureVector exons;
        transcriptTree->getMatching(exons, [](const GxfFeature* f) {
                return f->fType == GxfFeature::EXON;
            });
        return exons;
    }

    /* build transcript exons PSL to query and mapping to target genome.
     * Return NULL if no mappings for whatever reason.*/
    static PslMapping* exonTransMap(const TransMap* genomeTransMap,
                                    const FeatureNode* transcriptTree) {
        const string& qName(transcriptTree->fFeature->getAttr("transcript_id")->fVal);
        GxfFeatureVector exons = getExons(transcriptTree);
        // get alignment of exons to srcGenome and to targetGenome
        PslMapping* exonsMapping = FeatureTransMap(genomeTransMap).mapFeatures(qName, exons);
        if (debug) {
            exonsMapping->dump(cerr, "Transcript Exons:", "    ");
        }
        if ((exonsMapping == NULL) or (not exonsMapping->haveMappings())) {
            delete exonsMapping;
            return NULL;
        } else {
            return exonsMapping;
        }
    }
    
    /* create transMap objects used to do two level mapping via exons. */
    static const TransMapVector makeViaExonsTransMap(const PslMapping* exonsMapping) {
        TransMapVector transMaps;
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fSrcPsl, true)); // swap map genomeA to exons
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fMappedPsls[0], false)); // exons to genomeB
        return transMaps;
    }

    /* create a new transcript record that covers the alignment */
    void mapTranscriptFeature(FeatureNode* transcriptNode) {
        struct psl* mappedPsl = fExonsMapping->fMappedPsls[0];
        // transcript for mapped PSLs
        FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping,
                                   string(mappedPsl->tName),
                                   mappedPsl->tStart, mappedPsl->tEnd,
                                   charToString(pslQStrand(mappedPsl)));
        transcriptNode->fNumMappings = fExonsMapping->fMappedPsls.size();
        // if any parts was unmapped, also need a copy of the original transcript
        if (not pslQueryFullyMapped(mappedPsl)) {
            FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping);
        }
    }
    
    /* recursive map features below transcript */
    void mapFeatures(FeatureNode* featureNode) {
        const AttrVal* idAttr = featureNode->fFeature->findAttr("ID");
        const string& nodeId = (idAttr != NULL) ? idAttr->fVal : "someFeature";
        PslMapping* pslMapping = fViaExonsFeatureTransMap->mapFeature(nodeId, featureNode->fFeature);
        FeatureMapper::map(featureNode, pslMapping, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
           mapFeatures(featureNode->fChildren[iChild]);
        }
        delete pslMapping;
    }
    
    /* do work of mapping features when we have transmap mapping alignments
     * and know something will map */
    void mapTranscriptFeatures(FeatureNode* transcriptTree) {
        mapTranscriptFeature(transcriptTree);
        for (int iChild = 0; iChild < transcriptTree->fChildren.size(); iChild++) {
            mapFeatures(transcriptTree->fChildren[iChild]);
        }
    }

    /* recursive process features that are unmapped */
    void processUnmappedFeatures(FeatureNode* featureNode) {
        FeatureMapper::map(featureNode, NULL, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
            processUnmappedFeatures(featureNode->fChildren[iChild]);
        }
    }

    public:
    /* constructor */
    TranscriptMapper(const TransMap* genomeTransMap,
                     FeatureNode* transcriptTree,
                     bool srcSeqInMapping):
        fSrcSeqInMapping(srcSeqInMapping),
        fExonsMapping(exonTransMap(genomeTransMap, transcriptTree)),
        fViaExonsFeatureTransMap(NULL) {
        if (fExonsMapping != NULL) {
            fVaiExonsTransMaps = makeViaExonsTransMap(fExonsMapping);
            fViaExonsFeatureTransMap = new FeatureTransMap(fVaiExonsTransMaps);
        }
    }

    /* destructor */
    ~TranscriptMapper() {
        delete fExonsMapping;
        fVaiExonsTransMaps.free();
        delete fViaExonsFeatureTransMap;
    }
    
    /*
     * map one transcript's annotations.  Fill in transcriptTree
     */
    void mapTranscript(FeatureNode* transcriptTree) {
        assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
        if (fViaExonsFeatureTransMap != NULL) {
            mapTranscriptFeatures(transcriptTree);
        } else {
            processUnmappedFeatures(transcriptTree);
        }
        FeatureMapper::updateIds(transcriptTree);
        transcriptTree->setNumMappingsAttr();
    }
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
void GeneMapper::processTranscript(FeatureNode* transcriptTree) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcriptTree,
                                      isSrcSeqInMapping(transcriptTree));
    transcriptMapper.mapTranscript(transcriptTree);
}

/* process all transcripts of gene. */
void GeneMapper::processTranscripts(FeatureNode* geneTree) const {
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        processTranscript(transcriptTree);
    }
}

/* Recursively made features fulled unmapped, removing mapped and partially
 * mapped entries.  Used when we discover conflicts at the gene level. */
void GeneMapper::forceToUnmapped(FeatureNode* featureNode,
                                 RemapStatus remapStatus) const {
    FeatureMapper::forceToUnmapped(featureNode, remapStatus);
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
    FeatureMapper::mapBounding(geneTree, isSrcSeqInMapping(geneTree), seqid, start, end, strand);
}

/* If there are any unmapped transcripts of a gene, add a gene
 * record for the bounds */
void GeneMapper::buildUnmappedGeneFeature(FeatureNode* geneTree,
                                          bool srcSeqInMapping) const {
    FeatureMapper::mapBounding(geneTree, isSrcSeqInMapping(geneTree));
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
                              ostream& mappedGxfFh) {
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
                                ostream& unmappedGxfFh) {
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
void GeneMapper::outputInfoHeader(ostream& mappingInfoFh) {
    for (int i = 0; mappingInfoHeaders[i] != NULL; i++) {
        if (i > 0) {
            mappingInfoFh << "\t";
        }
        mappingInfoFh << mappingInfoHeaders[i];
    }
    mappingInfoFh << endl;
}

/* output info on one bounding feature */
void GeneMapper::outputFeatureInfo(FeatureNode* featureNode,
                                   ostream& mappingInfoFh) {
    assert(featureNode->fMappedFeatures.size() <= 1);  // only handles bounding features
    mappingInfoFh << featureNode->fFeature->getTypeId() << "\t"
                  << featureNode->fFeature->fType << "\t"
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
                  << featureNode->fNumMappings << endl;
}

/*
 * Output information about gene mapping
 */
void GeneMapper::outputInfo(FeatureNode* geneNode,
                            ostream& mappingInfoFh) {
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
                             ostream& mappingInfoFh) {
    FeatureNode* geneTree = GeneTree::geneTreeFactory(gxfParser, geneFeature);
    processTranscripts(geneTree);
    buildGeneFeature(geneTree);
    // handle gene level conflicted */
    if (hasMixedMappedSeqStrand(geneTree)) {
        forceToUnmapped(geneTree, REMAP_STATUS_GENE_CONFLICT);
    } else if (hasExcessiveSizeChange(geneTree)) {
        forceToUnmapped(geneTree, REMAP_STATUS_GENE_SIZE_CHANGE);
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
                               ostream& mappingInfoFh) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        processGene(gxfParser, dynamic_cast<GxfFeature*>(gxfRecord), mappedGxfFh, unmappedGxfFh, mappingInfoFh);
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
                        ostream& mappingInfoFh) {
    outputInfoHeader(mappingInfoFh);
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        processRecord(gxfParser, gxfRecord, mappedGxfFh, unmappedGxfFh, mappingInfoFh);
    }
}

