#include "geneMapper.hh"
#include "featureMapping.hh"
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
 * class to map a single transcript and subfeatures
 */
class TranscriptMapper {
    private:
    const GxfFeatureNode* fTranscriptTree;  // not owned
    const string& fQName;
    const TransMap* fGenomeTransMap;  // not owned
    bool fSrcSeqInMapping;  // is the source sequence in the genome mappings?
    FeatureMappingSet* fExonMappings;
    
    /* get exon features */
    GxfFeatureVector getExons() const {
        GxfFeatureVector exons;
        fTranscriptTree->getMatching(exons, [](const GxfFeature* f) {return f->fType == GxfFeature::EXON;});
        return exons;
    }

    /* map all exons of a transcript, store in fExonMappings */
    bool mapExons() {
        GxfFeatureVector exons = getExons();
        FeatureTransMap exonTransMap(fGenomeTransMap);
        PslMapping* pslMapping = exonTransMap.mapFeatures(fQName, exons);
        fExonMappings = FeatureMapper::mapFeatures(exons, pslMapping, fSrcSeqInMapping);
        return fSrcSeqInMapping && pslMapping->haveMappings();
    }

    /* process all unmapped features besides exons */
    
    public:
    /* constructor */
    TranscriptMapper(const GxfFeatureNode* transcriptTree,
                     const TransMap* genomeTransMap):
        fTranscriptTree(transcriptTree),
        fQName(transcriptTree->fFeature->getAttr("transcript_id")->fVal),
        fGenomeTransMap(genomeTransMap),
        fSrcSeqInMapping(genomeTransMap->haveTargetSeq(transcriptTree->fFeature->fSeqid)),
        fExonMappings(NULL) {
    }

    /* destructor */
    ~TranscriptMapper() {
        delete fExonMappings;
    }
    
    /*
     * map one transcript's annotations.
     */
    void mapTranscript() {
        if (mapExons()) {
        } else {
        }
    }
};

/* process one transcript */
void GeneMapper::processTranscript(const GxfFeatureNode* transcriptTree) const {
    TranscriptMapper transcriptMapper(transcriptTree, fGenomeTransMap);
    transcriptMapper.mapTranscript();
}

/*
 * map and output one gene's annotations
 */
void GeneMapper::processGene(GxfParser *gxfParser,
                             const GxfFeature* geneFeature,
                             ostream& outFh) const {
    GxfFeatureTree* geneTree = new GxfFeatureTree(gxfParser, geneFeature);
    for (size_t i = 0; i < geneTree->fGene->fChildren.size(); i++) {
        const GxfFeatureNode* transcriptTree = geneTree->fGene->fChildren[i];
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
