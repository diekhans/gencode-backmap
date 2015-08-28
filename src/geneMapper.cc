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

// FIXME: tmp
#define debug 0

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
    const FeatureTransMap* fViaExonsTransMap;   // two-level transmap, NULL if can't map (owned)
    
    /* get exon features */
    static GxfFeatureVector getExons(const FeatureNode* transcriptTree) {
        GxfFeatureVector exons;
        transcriptTree->getMatching(exons, [](const GxfFeature* f) {
                return f->fType == GxfFeature::EXON;
            });
        return exons;
    }

    /*
     * build exon PSL to query and mapping to target genome.  Return NULL if
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
    
    /* map transMap of exons of a transcript to the genome. */
    static const FeatureTransMap* makeViaExonsTransMap(const PslMapping* exonsMapping) {
        if (debug) {
            cerr << "exonsSrc:\t" << pslToString(exonsMapping->fSrcPsl) << endl;
            for (int i = 0; i < exonsMapping->fMappedPsls.size(); i++) {
                cerr << "exonsMap[" << i <<"]\t" << pslToString(exonsMapping->fMappedPsls[i]) << endl;
            }
        }
        TransMapVector transMaps;
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fSrcPsl, true)); // swap map genomeA to exons
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fMappedPsls[0], false)); // exons to genomeB
        return new FeatureTransMap(transMaps);
    }

    /* create a new transcript record that covers the alignment */
    void mapTranscriptFeature(FeatureNode* transcriptNode) {
        struct psl* mappedPsl = fExonsMapping->fMappedPsls[0];
        // transcript for mapped PSLs
        FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping,
                                   string(mappedPsl->tName),
                                   mappedPsl->tStart, mappedPsl->tEnd,
                                   charToString(pslQStrand(mappedPsl)));
        // if any was unmapped, also need a copy of the original transcript
        if (not pslQueryFullyMapped(mappedPsl)) {
            FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping);
        }

        
    }
    
    /* recursive map features below transcript */
    void mapFeatures(FeatureNode* featureNode) {
        PslMapping* pslMapping = fViaExonsTransMap->mapFeature("someFeature", featureNode->fFeature);
        if (debug && pslMapping != NULL) {
            cerr << "src\t" << pslToString(pslMapping->fSrcPsl) << endl;
            if (pslMapping->haveMappings()) {
                cerr << "mapped\t" << pslToString(pslMapping->fMappedPsls[0]) << endl;
            }
        }
        FeatureMapper::map(featureNode, pslMapping, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
           mapFeatures(featureNode->fChildren[iChild]);
        }
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
        fViaExonsTransMap(NULL) {
        if (fExonsMapping != NULL) {
            fViaExonsTransMap = makeViaExonsTransMap(fExonsMapping);
        }
    }

    /* destructor */
    ~TranscriptMapper() {
        delete fExonsMapping;
        delete fViaExonsTransMap;
    }
    
    /*
     * map one transcript's annotations.  Fill in transcriptTree
     */
    void mapTranscript(FeatureNode* transcriptTree) {
        assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
        if (fViaExonsTransMap != NULL) {
            mapTranscriptFeatures(transcriptTree);
        } else {
            processUnmappedFeatures(transcriptTree);
        }
        FeatureMapper::updateIds(transcriptTree);
    }
};

/* is the source sequence for a feature in the mapping at all? */
bool GeneMapper::isSrcSeqInMapping(FeatureNode* featureNode) const {
    return fGenomeTransMap->haveTargetSeq(featureNode->fFeature->fSeqid);
}

/* process one transcript */
void GeneMapper::processTranscript(FeatureNode* transcriptTree) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcriptTree,
                                      isSrcSeqInMapping(transcriptTree));
    transcriptMapper.mapTranscript(transcriptTree);
}

/* process all transcripts of gene. */
void GeneMapper::processTranscripts(FeatureNode* geneNode) const {
    for (size_t i = 0; i < geneNode->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneNode->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        processTranscript(transcriptTree);
    }
}

/* are there any mapped transcripts? */
bool GeneMapper::haveMappedTranscripts(FeatureNode* geneNode) const {
    return (geneNode->fChildren.size() > 0) and (geneNode->fChildren[0]->fMappedFeatures.size() > 0);
}

/* are there any unmapped transcripts? */
bool GeneMapper::haveUnmappedTranscripts(FeatureNode* geneNode) const {
    return (geneNode->fChildren.size() > 0) and (geneNode->fChildren[0]->fUnmappedFeatures.size() > 0);
}

/* If there are any mapped transcripts of a gene, add a gene
 * record for the bounds */
void GeneMapper::buildMappedGeneFeature(FeatureNode* geneNode,
                                        bool srcSeqInMapping) const {
    // first transcripts bounds
    FeatureNode* transcriptTree = geneNode->fChildren[0];
    assert(transcriptTree->fMappedFeatures.size() == 1);
    GxfFeature* feature = transcriptTree->fMappedFeatures[0];
    const string& seqid = feature->fSeqid;
    const string& strand = feature->fStrand;
    int start = feature->fStart, end = feature->fEnd;
    // other transcripts bounds
    for (int i = 1; i < geneNode->fChildren.size(); i++) {
        transcriptTree = geneNode->fChildren[i];
        assert(transcriptTree->fMappedFeatures.size() == 1);
        feature = transcriptTree->fMappedFeatures[0];
        if (feature->fSeqid != seqid) {
            throw logic_error("gene has transcripts on different sequences: " + feature->toString());
        }
        if (feature->fStrand != strand) {
            throw logic_error("gene has transcripts on different strands: " + feature->toString());
        }
        start = min(start, feature->fStart);
        end = max(end, feature->fEnd);
    }
    FeatureMapper::mapBounding(geneNode, isSrcSeqInMapping(transcriptTree), seqid, start, end, strand);
}

/* If there are any unmapped transcripts of a gene, add a gene
 * record for the bounds */
void GeneMapper::buildUnmappedGeneFeature(FeatureNode* geneNode,
                                          bool srcSeqInMapping) const {
    FeatureMapper::mapBounding(geneNode, isSrcSeqInMapping(geneNode));
}

/* output GFF3 mapped ##sequence-region if not already written */
void GeneMapper::outputMappedSeqRegionIfNeed(const GxfFeature* feature,
                                             ostream& mappedGxfFh) {
    if ((feature->getFormat() == GFF3_FORMAT)
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

/*
 * map and output one gene's annotations
 */
void GeneMapper::processGene(GxfParser *gxfParser,
                             GxfFeature* geneFeature,
                             ostream& mappedGxfFh,
                             ostream& unmappedGxfFh) {
    FeatureTree* geneTree = new FeatureTree(gxfParser, geneFeature);
    processTranscripts(geneTree->fGene);
    bool srcSeqInMapping = fGenomeTransMap->haveTargetSeq(geneTree->fGene->fFeature->fSeqid);
    if (haveMappedTranscripts(geneTree->fGene)) {
        buildMappedGeneFeature(geneTree->fGene, srcSeqInMapping);
        outputMappedSeqRegionIfNeed(geneTree->fGene->fFeature, mappedGxfFh);
        outputMapped(geneTree->fGene, mappedGxfFh);
    }
    if (haveUnmappedTranscripts(geneTree->fGene)) {
        buildUnmappedGeneFeature(geneTree->fGene, srcSeqInMapping);
        outputUnmapped(geneTree->fGene, unmappedGxfFh);
    }
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

