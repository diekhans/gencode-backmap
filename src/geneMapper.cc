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

/* fraction of gene expansion that causes a rejection */
const float geneExpansionThreshold = 0.20;  


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

    /* build exon PSL to query and mapping to target genome.  Return NULL if
     * no mappings for whatever reason.*/
    static PslMapping* exonTransMap(const TransMap* genomeTransMap,
                                    const FeatureNode* transcriptTree) {
        const string& qName(transcriptTree->fFeature->getAttr("transcript_id")->fVal);
        GxfFeatureVector exons = getExons(transcriptTree);
        // get alignment of exons to srcGenome and to targetGenome
        PslMapping* exonsMapping = FeatureTransMap(genomeTransMap).mapFeatures(qName, exons);
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
                                   charToString(pslQStrand(mappedPsl)),
                                   fExonsMapping->fMappedPsls.size());
        // if any was unmapped, also need a copy of the original transcript
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
    }
};

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(const GxfFeature* feature) const {
    return fGenomeTransMap->haveTargetSeq(feature->fSeqid);
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
void GeneMapper::forceTranscriptToUnmapped(FeatureNode* featureNode,
                                           RemapStatus remapStatus) const {
    FeatureMapper::forceToUnmapped(featureNode, remapStatus);
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        forceTranscriptToUnmapped(featureNode->fChildren[i], remapStatus);
    }
}

/* Force all transcript features to fulled unmapped.  Used when we discover
 * conflicts at the gene level. */
void GeneMapper::forceTranscriptsToUnmapped(FeatureNode* geneTree,
                                            RemapStatus remapStatus) const {
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        forceTranscriptToUnmapped(geneTree->fChildren[i], remapStatus);
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

/* check if gene transcripts cause the gene to expand it's length
 * beyond a threshold*/
bool GeneMapper::hasExcessiveExpansion(FeatureNode* geneTree) const {
    int srcGeneLength = geneTree->fFeature->size();
    int mappedGeneLength = calcMappedGeneLength(geneTree);
    if (mappedGeneLength > srcGeneLength) {
        return (float(mappedGeneLength-srcGeneLength) / float(srcGeneLength)) > geneExpansionThreshold;
    } else {
        return false;
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
void GeneMapper::buildGeneFeatures(FeatureNode* geneTree) const {
    bool srcSeqInMapping = fGenomeTransMap->haveTargetSeq(geneTree->fFeature->fSeqid);
    if (haveMappedTranscripts(geneTree)) {
        buildMappedGeneFeature(geneTree, srcSeqInMapping);
    }
    if (haveUnmappedTranscripts(geneTree)) {
        buildUnmappedGeneFeature(geneTree, srcSeqInMapping);
    }
}

/* output GFF3 mapped ##sequence-region if not already written */
void GeneMapper::outputMappedSeqRegionIfNeed(const GxfFeature* feature,
                                             ostream& mappedGxfFh) {
    if ((feature->getFormat() == GFF3_FORMAT)
        and isSrcSeqInMapping(feature)
        and (not checkRecordSeqRegionWritten(feature->fSeqid))) {
        mappedGxfFh << "##sequence-region " << feature->fSeqid << " 1 "
                    << fGenomeTransMap->getTargetSeqSize(feature->fSeqid) << endl;
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
void GeneMapper::output(FeatureNode* geneNode,
                        ostream& mappedGxfFh,
                        ostream& unmappedGxfFh) {
    if (haveMappedTranscripts(geneNode)) {
        outputMappedSeqRegionIfNeed(geneNode->fFeature, mappedGxfFh);
        outputMapped(geneNode, mappedGxfFh);
    }
    if (haveUnmappedTranscripts(geneNode)) {
        outputUnmapped(geneNode, unmappedGxfFh);
    }
}

/*
 * map and output one gene's annotations
 */
void GeneMapper::processGene(GxfParser *gxfParser,
                             GxfFeature* geneFeature,
                             ostream& mappedGxfFh,
                             ostream& unmappedGxfFh) {
    FeatureTree* geneTree = new FeatureTree(gxfParser, geneFeature);
    FeatureNode* geneNode = geneTree->fGene;
    processTranscripts(geneNode);
    // handle gene level conflicted */
    if (hasMixedMappedSeqStrand(geneNode)) {
        forceTranscriptsToUnmapped(geneNode, REMAP_STATUS_GENE_CONFLICT);
    }
    if (hasExcessiveExpansion(geneNode)) {
        forceTranscriptsToUnmapped(geneNode, REMAP_STATUS_GENE_EXPAND);
    }
    buildGeneFeatures(geneNode);
    geneNode->setRemapStatusAttr();
    output(geneNode, mappedGxfFh, unmappedGxfFh);
    delete geneTree;
}

/* process a record, this may consume additional feature records  */
void GeneMapper::processRecord(GxfParser *gxfParser,
                               GxfRecord* gxfRecord,
                               ostream& mappedGxfFh,
                               ostream& unmappedGxfFh) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        processGene(gxfParser, dynamic_cast<GxfFeature*>(gxfRecord), mappedGxfFh, unmappedGxfFh);
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
                        ostream& unmappedGxfFh) {
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        processRecord(gxfParser, gxfRecord, mappedGxfFh, unmappedGxfFh);
    }
}

