
#include "geneMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "featureTransMap.hh"
#include <fstream>
#include <iostream>

/* get exon features */
GxfFeatureVector GeneMapper::getExons(const GxfFeatureNode* transcriptTree) const {
    GxfFeatureVector exons;
    for (size_t i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            exons.push_back(transcriptTree->fChildren[i]->fFeature);
        }
    }
    return exons;
}

/* 
 * Map a transcript's exons, returning object that has mappings and scoring
 */
PslMapping* GeneMapper::mapTranscriptExons(const GxfFeatureNode* transcriptTree) const {
    const string& qName = transcriptTree->fFeature->getAttr("transcript_id")->fVal;
    GxfFeatureVector exons = getExons(transcriptTree);
    return fFeatureTransMap->mapFeatures(qName, exons);
}

/*
 * map and output one transcripts annotations.
 */
void GeneMapper::processTranscript(const GxfFeatureNode* transcriptTree,
                                   ostream& outFh) const {
    PslMapping *pslMapping = mapTranscriptExons(transcriptTree);
    transcriptTree->write(outFh);
    outFh << "srcPsl: " << pslToString(pslMapping->fSrcPsl) << endl;
    outFh <<"destPsls: score: " << pslMapping->fScore
          << " srcBlks: " << pslMapping->fSrcPsl->blockCount
          << " destBlks: " << ((pslMapping->fMappedPsls.size() == 0) ? 0 : pslMapping->fMappedPsls[0]->blockCount)
          << endl;
    for (size_t i = 0; i < pslMapping->fMappedPsls.size(); i++) {
        outFh << "        " << pslToString(pslMapping->fMappedPsls[i]) << endl;
    }
    outFh << endl;
    delete pslMapping;
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
            throw invalid_argument("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        processTranscript(transcriptTree, outFh);
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
            outFh << gxfRecord->toString() << endl;
            delete gxfRecord;
        }
    }
}
