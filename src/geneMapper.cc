
#include "geneMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "pslOps.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include "Frame.hh"
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
     * create an exon feature for a full or partially mapped feature.
     */
    const GxfFeature* mkMappedFeature(const GxfFeature* exon,
                                      const PslCursor& srcPslCursor,
                                      const PslCursor& mappedPslCursor,
                                      int length,
                                      RemapStatus remapStatus) {
        int off = srcPslCursor.getTPosStrand('+') - (exon->fStart-1);
        Frame frame(Frame::fromPhaseStr(exon->fPhase).incr(off));

        // GFF3 genomic coordinates are always plus strand
        int mappedTStart = mappedPslCursor.getTPosStrand('+');
        int mappedTEnd = mappedTStart + length;

        return gxfFeatureFactory(exon->getFormat(), exon->fSeqid, exon->fSource, exon->fType,
                                 mappedTStart, mappedTEnd, exon->fScore,
                                 string(0, pslQStrand(srcPslCursor.getPsl())),
                                 frame.toPhaseStr(), exon->fAttrs, &remapStatusToAttrVal(remapStatus));
    }
    
    /* 
     * Create an exon feature for a full or partially unmapped feature.
     * The feature is in source coordinates.
     */
    const GxfFeature* mkUnmappedFeature(const GxfFeature* exon,
                                        const PslCursor& srcPslCursor,
                                        const PslCursor& mappedPslCursor,
                                        int length,
                                        RemapStatus remapStatus) {
        int off = srcPslCursor.getTPosStrand('+') - (exon->fStart-1);
        Frame frame(Frame::fromPhaseStr(exon->fPhase).incr(off));

        // GFF3 genomic coordinates are always plus strand
        int unmappedTStart = srcPslCursor.getTPosStrand('+');
        int unmappedTEnd = unmappedTStart + length;

        return gxfFeatureFactory(exon->getFormat(), exon->fSeqid, exon->fSource, exon->fType,
                                 unmappedTStart, unmappedTEnd, exon->fScore, exon->fStrand,
                                 frame.toPhaseStr(), exon->fAttrs, &remapStatusToAttrVal(remapStatus));
    }

    /* compute length of mapped region of feature */
    static int calcMappedLength(const PslCursor& srcPslCursor,
                                const PslCursor& mappedPslCursor) {
        // length is the minimum left in either block
        return min(srcPslCursor.getBlockLeft(), mappedPslCursor.getBlockLeft());
    }

    /* compute length of unmapped internal region */
    static int calcInternalUnmappedLength(const PslCursor& srcPslCursor,
                                          const PslCursor& mappedPslCursor) {
        // length is minimum of different between starts in feature and how much is left in the feature
        return min(mappedPslCursor.getQPos()-srcPslCursor.getQPos(), srcPslCursor.getBlockLeft());
    }

    /* compute length of unmapped terminal region */
    static int calcTerminalUnmappedLength(const PslCursor& srcPslCursor,
                                          const PslCursor& mappedPslCursor) {
        //  length is what is left over in this src block
        return srcPslCursor.getBlockLeft();
    }

   
    /*
     * Determine status of one part of a feature.  Return true if
     * the part was mapped, false if it wasn't. Cursors are updated
     */
    boolean featureMappingStatusPart(const GxfFeature* feature,
                                     PslCursor& srcPslCursor,
                                     PslCursor& mappedPslCursor) const {
        assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
        if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
            // deleted region
            int length = calcInternalUnmappedLength(srcPslCursor, mappedPslCursor);
            srcPslCursor = srcPslCursor.advance(length);
            return false;
        } else {
            // mapped region
            int length = calcMappedLength(srcPslCursor, mappedPslCursor);
            srcPslCursor = srcPslCursor.advance(length);
            mappedPslCursor = mappedPslCursor.advance(length);
            return true;
        }
    }

    /* Determine mapping status of an feature.  This parallels the logic in
     * mapExon, we look at the mapping twice so we have the status to
     * create the new features as we go.
     */
    RemapStatus featureMappingStatus(const GxfFeature* feature,
                                     const PslCursor& srcPslCursorIn,
                                     const PslCursor& mappedPslCursorIn) const {
        PslCursor srcPslCursor(srcPslCursorIn);
        PslCursor mappedPslCursor(mappedPslCursorIn);
        int mapCnt = 0, deleteCnt = 0;

        int srcPslFeatureQEnd = srcPslCursor.getQBlockEnd();
        while ((srcPslCursor.getQPos() < srcPslFeatureQEnd) and (not mappedPslCursor.atEnd())) {
            if (featureMappingStatusPart(feature, srcPslCursor, mappedPslCursor)) {
                mapCnt++;
            } else {
                deleteCnt++;
            }
        }
        if (srcPslCursor.getQPos() < srcPslFeatureQEnd) {
            deleteCnt++;
        }
        if (deleteCnt == 0) {
            // fully mapped
            return (mapCnt == 1) ? REMAP_STATUS_FULL_CONTIG : REMAP_STATUS_FULL_FRAGMENT;
        } else if (mapCnt == 0) {
            // feature not mapped
            return REMAP_STATUS_DELETED;
        } else {
            // partially mapped
            return (mapCnt == 1) ? REMAP_STATUS_PARTIAL_CONTIG : REMAP_STATUS_PARTIAL_FRAGMENT;
        }
    }

    /*
     * Map one part of an exon feature.  Cursors are updated
     */
    void mapExonPart(const GxfFeature* exon,
                     PslCursor& srcPslCursor,
                     PslCursor& mappedPslCursor,
                     RemapStatus remapStatus) {
        assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
        if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
            // deleted region
            int length = calcInternalUnmappedLength(srcPslCursor, mappedPslCursor);
            delete mkUnmappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus);
            srcPslCursor = srcPslCursor.advance(length);
        } else {
            // mapped region
            int length = calcMappedLength(srcPslCursor, mappedPslCursor);
            delete mkMappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus);
            srcPslCursor = srcPslCursor.advance(length);
            mappedPslCursor = mappedPslCursor.advance(length);
        }
    }

    /*
     * Map an exon feature.  Cursors are updated and srcPslCursor will point to the
     * next feature.  This parallels logic in featureMappingStatus.
     */
    void mapExon(const GxfFeature* exon,
                 PslCursor& srcPslCursor,
                 PslCursor& mappedPslCursor,
                 RemapStatus remapStatus) {
        assert((exon->fEnd-exon->fStart)+1 == srcPslCursor.getBlockLeft());

        // note that source blocks can be merged in mapped block, but we don't merge
        // features.
        int srcPslExonQEnd = srcPslCursor.getQBlockEnd();
        while ((srcPslCursor.getQPos() < srcPslExonQEnd) && (not mappedPslCursor.atEnd())) {
            mapExonPart(exon, srcPslCursor, mappedPslCursor, remapStatus);
        }
        if (srcPslCursor.getQPos() < srcPslExonQEnd) {
            // unmapped at the end of feature
            int length = calcTerminalUnmappedLength(srcPslCursor, mappedPslCursor);
            delete mkUnmappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus);
            srcPslCursor = srcPslCursor.advance(length);
        }
        assert(srcPslCursor.getQPos() == srcPslExonQEnd);
    }

    /*
     * Map exons of a transcript that has at least something mapped.
     */
    void mapExons(const GxfFeatureVector& exons) {
        PslCursor srcPslCursor(fPslMapping->fSrcPsl);
        PslCursor mappedPslCursor(fPslMapping->fMappedPsls[0]);
        for (int iExon = 0; iExon < exons.size(); iExon++) {
            RemapStatus remapStatus = featureMappingStatus(exons[iExon], srcPslCursor, mappedPslCursor);
            mapExon(exons[iExon], srcPslCursor, mappedPslCursor, remapStatus);
        }
    }

    /* map all exons of a transcript mapping failure */
    bool mapTranscriptExons(const FeatureTransMap* featureTransMap) {
        GxfFeatureVector exons = getExons(fTranscriptTree);
        fPslMapping = featureTransMap->mapFeatures(fQName, exons);
        if (fPslMapping == NULL) {
            return false;
        } else if (fPslMapping->fMappedPsls.size() == 0) {
            return false;
        } else {
            mapExons(exons);
            return true;
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
