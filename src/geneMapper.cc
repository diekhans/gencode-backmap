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
    static GxfFeatureVector getExons(const GxfFeatureNode* transcriptTree) {
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
                                    const GxfFeatureNode* transcriptTree) {
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
        TransMapVector transMaps;
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fSrcPsl, true)); // swap map genomeA to exons
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fMappedPsls[0], false)); // exons to genomeB
        return new FeatureTransMap(transMaps);
    }

    /* create a new transcript record that covers the alignment */
    void mapTranscriptFeature(GxfFeatureNode* transcriptNode) {
        struct psl* mappedPsl = fExonsMapping->fMappedPsls[0];
        FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping,
                                   string(mappedPsl->tName),
                                   mappedPsl->tStart+1, mappedPsl->tEnd,
                                   pslQStrand(mappedPsl));
    }
    
    /* recursive map features below transcript */
    void mapFeatures(GxfFeatureNode* featureNode) {
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
            processUnmappedFeatures(featureNode->fChildren[iChild]);
        }
    }
    
    /* do work of mapping features when we have transmap mapping alignments
     * and know something will map */
    void mapTranscriptFeatures(GxfFeatureNode* transcriptTree) {
        mapTranscriptFeature(transcriptTree);
        for (int iChild = 0; iChild < transcriptTree->fChildren.size(); iChild++) {
            mapFeatures(transcriptTree->fChildren[iChild]);
        }
    }
    
    /* recursive process features that are unmapped */
    void processUnmappedFeatures(GxfFeatureNode* featureNode) {
        FeatureMapper::map(featureNode, NULL, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
            processUnmappedFeatures(featureNode->fChildren[iChild]);
        }
    }

    public:
    /* constructor */
    TranscriptMapper(const TransMap* genomeTransMap,
                     GxfFeatureNode* transcriptTree):
        fSrcSeqInMapping(genomeTransMap->haveTargetSeq(transcriptTree->fFeature->fSeqid)),
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
    void mapTranscript(GxfFeatureNode* transcriptTree) {
        assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
        if (fViaExonsTransMap != NULL) {
            mapTranscriptFeatures(transcriptTree);
        } else {
            processUnmappedFeatures(transcriptTree);
        }
    }
};

/* process one transcript */
void GeneMapper::processTranscript(GxfFeatureNode* transcriptTree) const {
    TranscriptMapper transcriptMapper(fGenomeTransMap, transcriptTree);
    transcriptMapper.mapTranscript(transcriptTree);
    transcriptTree->dump(cerr);
}

/*
 * map and output one gene's annotations
 */
void GeneMapper::processGene(GxfParser *gxfParser,
                             const GxfFeature* geneFeature,
                             ostream& outFh) const {
    GxfFeatureTree* geneTree = new GxfFeatureTree(gxfParser, geneFeature);
    for (size_t i = 0; i < geneTree->fGene->fChildren.size(); i++) {
        GxfFeatureNode* transcriptTree = geneTree->fGene->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        processTranscript(transcriptTree);
    }
    delete geneTree;
}
                             

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfParser *gxfParser,
                        ostream& outFh) const {
    const GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (instanceOf(gxfRecord, GxfFeature)) {
            processGene(gxfParser, dynamic_cast<const GxfFeature*>(gxfRecord), outFh);
        } else {
            delete gxfRecord;
        }
    }
}
