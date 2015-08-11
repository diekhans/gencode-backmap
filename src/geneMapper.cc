
#include "geneMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "pslOps.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include <fstream>
#include <iostream>

/**
 * class to map a single transcript and subfeatures
 */
class TranscriptMapper {
    private:
    const GxfFeatureNode* fTranscriptTree;
    const string fQName;
    PslMapping *fPslMapping;

    
    /* get exon features */
    GxfFeatureVector getExons(const GxfFeatureNode* transcriptTree) const {
        GxfFeatureVector exons;
        for (size_t i = 0; i < transcriptTree->fChildren.size(); i++) {
            if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
                exons.push_back(transcriptTree->fChildren[i]->fFeature);
            }
        }
        return exons;
    }

    /* 
     * create an exon feature for a full mapped
     */
    const GxfFeature* mkFullMappedFeature(const GxfFeature* exon,
                                          const PslCursor& srcPslCursor,
                                          const PslCursor& mappedPslCursor,
                                          int amount) {
        assert(srcPslCursor.getTPosStrand('+') == exon->fStart-1);
        assert(srcPslCursor.getTBlockEndStrand('+') == exon->fEnd);
        assert(srcPslCursor.getTPos()+amount == exon->fEnd);

        // coordinates are always plus strand
        int mappedTStart = mappedPslCursor.getTPosStrand('+');
        int mappedTEnd = mappedTStart + amount;

        return gxfFeatureFactory(exon->getFormat(), exon->fSeqid, exon->fSource, exon->fType,
                                 mappedTStart, mappedTEnd, exon->fScore,
                                 string(0, pslQStrand(srcPslCursor.getPsl())), ".", exon->fAttrs);
    }
    
    /*
     * Map one part of an exon feature.  Cursors are updated
     */
    void mapExonPart(const GxfFeature* exon,
                     PslCursor& srcPslCursor,
                     PslCursor& mappedPslCursor,
                     ostream& outFh) {
        assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
        outFh << "  exonPart: " << srcPslCursor.toString() << " = " << mappedPslCursor.toString() << endl;
        if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
            // deleted region
            int amount = min(mappedPslCursor.getQPos()-srcPslCursor.getQPos(), srcPslCursor.getBlockLeft());
            outFh << "    deleted: " << amount << endl;
            srcPslCursor = srcPslCursor.advance(amount);
        } else {
            // mapped region
            int amount = min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
            if (amount == exon->size()) {
                // whole 
                delete mkFullMappedFeature(exon, srcPslCursor, mappedPslCursor, amount);
            }
            outFh << "    overlap: " << amount << endl;
            srcPslCursor = srcPslCursor.advance(amount);
            mappedPslCursor = mappedPslCursor.advance(amount);
        }
    }

    /*
     * Map an exon feature.  Cursors are updated
     */
    void mapExon(const GxfFeature* exon,
                 PslCursor& srcPslCursor,
                 PslCursor& mappedPslCursor) {
        assert((exon->fEnd-exon->fStart)+1 == srcPslCursor.getBlockLeft());
        cerr << "mappingExon: "  << srcPslCursor.toString() << " = " << mappedPslCursor.toString() << "\t" << exon->toString() << endl;

        // note that source blocks can be merged in mapped block, but we don't merge
        // features
        int srcPslExonQEnd = srcPslCursor.getQBlockEnd();
        while ((srcPslCursor.getQPos() < srcPslExonQEnd) && (not mappedPslCursor.atEnd())) {
            mapExonPart(exon, srcPslCursor, mappedPslCursor, cerr);
        }
        // FIXME: handled un mapped at end
        if (srcPslCursor.getQPos() < srcPslExonQEnd) {
            cerr << "   exonLeftover: " << srcPslCursor.toString() << " = " << mappedPslCursor.toString() << endl;
        }
    }

    /*
     * Map exons of a transcript that has at least something mapped.
     */
    void mapExons(const GxfFeatureVector& exons,
                  PslMapping *pslMapping) {
        PslCursor srcPslCursor(pslMapping->fSrcPsl);
        PslCursor mappedPslCursor(pslMapping->fMappedPsls[0]);
        for (int iExon = 0; iExon < exons.size(); iExon++) {
            mapExon(exons[iExon], srcPslCursor, mappedPslCursor);
        }
    }

    /* map all exons of a transcript mapping failure */
    bool mapTranscriptExons(const FeatureTransMap* featureTransMap) {
        GxfFeatureVector exons = getExons(fTranscriptTree);
        fPslMapping = featureTransMap->mapFeatures(fQName, exons);
        if (fPslMapping->fMappedPsls.size() > 0) {
            mapExons(exons, fPslMapping);
            return true;
        } else {
            // outFh << "not mapped: " << qName << endl;
            return false;
        }
    }

    public:
    /* constructor */
    TranscriptMapper(const GxfFeatureNode* transcriptTree):
        fTranscriptTree(transcriptTree),
        fQName(transcriptTree->fFeature->getAttr("transcript_id")->fVal),
        fPslMapping(NULL) {
    }

    /* destructor */
    ~TranscriptMapper() {
        delete fPslMapping;
    }
    
    /*
     * map and output one transcripts annotations.
     */
    const GxfFeatureNode* mapTranscriptFeatures(const FeatureTransMap* featureTransMap) {
        mapTranscriptExons(featureTransMap);
        return NULL; // FIXME:
    }
};

/* process one transcript */
void GeneMapper::processTranscript(const GxfFeatureNode* transcriptTree) const {
    TranscriptMapper transcriptMapper(transcriptTree);
    transcriptMapper.mapTranscriptFeatures(fFeatureTransMap);
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
            // FIXME: outFh << gxfRecord->toString() << endl;
            delete gxfRecord;
        }
    }
}
