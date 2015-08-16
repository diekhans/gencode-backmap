
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

    GxfFeatureVector fMappedExons;
    GxfFeatureVector fUnmappedExons;


    /* record mapped exon */
    void recordMappedExon(const GxfFeature* feature) {
        fMappedExons.push_back(feature);
        cerr << "mapped\t" << feature->toString() << endl;
    }
    
    /* record unmapped exon */
    void recordUnmappedExon(const GxfFeature* feature) {
        fUnmappedExons.push_back(feature);
        cerr << "unmapped\t" << feature->toString() << endl;
    }
    
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

    /* update an id attribute if we have one, and make unique if necessary */
    static void updateIdAttr(AttrVals& attrs,
                             RemapStatus remapStatus,
                             int partIdx) {
        const AttrVal* idAttr = attrs.findAttr("ID");
        if ((idAttr != NULL) && ((remapStatus == REMAP_STATUS_FULL_FRAGMENT) or (remapStatus == REMAP_STATUS_PARTIAL_FRAGMENT))) {
            attrs.update(AttrVal(idAttr->fName, idAttr->fVal + "_" + toString(partIdx)));
        }
    }
    
    /* 
     * create an feature for a full or partially mapped feature.
     * partIdx is used to make ID unique if split.
     */
    const GxfFeature* mkMappedFeature(const GxfFeature* feature,
                                      const PslCursor& srcPslCursor,
                                      const PslCursor& mappedPslCursor,
                                      int length,
                                      RemapStatus remapStatus,
                                      int partIdx) {
        int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
        Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

        // GxF genomic coordinates are always plus strand
        int mappedTStart = mappedPslCursor.getTPosStrand('+');
        int mappedTEnd = mappedTStart + length;

        AttrVals attrs(feature->fAttrs);
        attrs.add(remapStatusToAttrVal(remapStatus));
        updateIdAttr(attrs, remapStatus, partIdx);
        
        return gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                                 mappedTStart, mappedTEnd, feature->fScore,
                                 string(0, pslQStrand(srcPslCursor.getPsl())),
                                 frame.toPhaseStr(), attrs);
    }
    
    /* 
     * Create an feature for a full or partially unmapped feature.
     * The created feature is in source coordinates.
     * partIdx is used to make ID unique if split.
     */
    const GxfFeature* mkUnmappedFeature(const GxfFeature* feature,
                                        const PslCursor& srcPslCursor,
                                        const PslCursor& mappedPslCursor,
                                        int length,
                                        RemapStatus remapStatus,
                                        int partIdx) {
        int off = srcPslCursor.getTPosStrand('+') - (feature->fStart-1);
        Frame frame(Frame::fromPhaseStr(feature->fPhase).incr(off));

        // GxF genomic coordinates are always plus strand
        int unmappedTStart = srcPslCursor.getTPosStrand('+');
        int unmappedTEnd = unmappedTStart + length;

        AttrVals attrs(feature->fAttrs);
        attrs.add(remapStatusToAttrVal(remapStatus));
        updateIdAttr(attrs, remapStatus, partIdx);
        
        return gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                                 unmappedTStart, unmappedTEnd, feature->fScore, feature->fStrand,
                                 frame.toPhaseStr(), attrs);
    }

    /* generate for unmapped feature for full unmapped transcript */
    const GxfFeature* mkUnmappedTranscriptFeature(const GxfFeature* feature,
                                                  RemapStatus remapStatus) {
        AttrVals attrs(feature->fAttrs);
        attrs.add(remapStatusToAttrVal(remapStatus));
        return gxfFeatureFactory(feature->getFormat(), feature->fSeqid, feature->fSource, feature->fType,
                                 feature->fStart, feature->fEnd, feature->fScore, feature->fStrand,
                                 feature->fPhase, attrs);
        
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
                     RemapStatus remapStatus,
                     int partIdx) {
        assert(srcPslCursor.getQPos() <= mappedPslCursor.getQPos());
        if (srcPslCursor.getQPos() < mappedPslCursor.getQPos()) {
            // deleted region
            int length = calcInternalUnmappedLength(srcPslCursor, mappedPslCursor);
            recordUnmappedExon(mkUnmappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus, partIdx));
            srcPslCursor = srcPslCursor.advance(length);
        } else {
            // mapped region
            int length = calcMappedLength(srcPslCursor, mappedPslCursor);
            recordMappedExon(mkMappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus, partIdx));
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
        int partIdx = 0;
        while ((srcPslCursor.getQPos() < srcPslExonQEnd) && (not mappedPslCursor.atEnd())) {
            mapExonPart(exon, srcPslCursor, mappedPslCursor, remapStatus, partIdx++);
        }
        if (srcPslCursor.getQPos() < srcPslExonQEnd) {
            // unmapped at the end of feature
            int length = calcTerminalUnmappedLength(srcPslCursor, mappedPslCursor);
            recordUnmappedExon(mkUnmappedFeature(exon, srcPslCursor, mappedPslCursor, length, remapStatus, partIdx++));
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

    /* handle none of the transcript was mapped */
    void processUnmappedTranscript(const GxfFeatureVector& exons) {
        // was there no mapping of entire sequence vs just not of this transcript?
        RemapStatus remapStatus = (fPslMapping == NULL) ? REMAP_STATUS_NO_SEQ_MAP : REMAP_STATUS_DELETED;

        for (int iExon = 0; iExon < exons.size(); iExon++) {
            recordUnmappedExon(mkUnmappedTranscriptFeature(exons[iExon], remapStatus));
        }
    }
    
    /* map all exons of a transcript mapping failure */
    bool mapTranscriptExons(const FeatureTransMap* featureTransMap) {
        GxfFeatureVector exons = getExons(fTranscriptTree);
        fPslMapping = featureTransMap->mapFeatures(fQName, exons);
        if ((fPslMapping == NULL) or (fPslMapping->fMappedPsls.size() == 0)) {
            processUnmappedTranscript(exons);
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
        fPslMapping(NULL),
        fMappedExons(true),
        fUnmappedExons(true) {
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
