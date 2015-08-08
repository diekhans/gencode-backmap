
#include "geneMapper.hh"
#include "jkinclude.hh"
#include "gxf.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include <fstream>
#include <iostream>

/*
 * Cursor into a PSL.  Tracks position in an alignment.
 */
class PslCursor {
    private:
    struct psl* fPsl;
    int fIBlk;  // index of current block, set to blockCount if reached end
    int fOff;   // offset in current block;

    /* check if in range */
    inline void assertInRange() const {
        assert(fIBlk < fPsl->blockCount);
        assert(fOff < fPsl->blockSizes[fIBlk]);
    }
    
    public:
    PslCursor(struct psl* psl, int iBlk=0, int off=0):
        fPsl(psl), fIBlk(iBlk), fOff(off) {
    }

    /* have we reached the end of the psl */
    bool atEnd() const {
        return fIBlk >= fPsl->blockCount;
    }

    /* accessors of current position */
    int getQPos() const {
        if (atEnd()) {
            return pslQEnd(fPsl, fPsl->blockCount-1); 
        } else {
            return pslQStart(fPsl, fIBlk)+fOff;
        }
    }
    int getTPos() const {
        if (atEnd()) {
            return pslTEnd(fPsl, fPsl->blockCount-1); 
        } else {
            return pslTStart(fPsl, fIBlk)+fOff;
        }
    }

    /* accessors for current block end */
    int getQBlockEnd() const {
        if (atEnd()) {
            return pslQEnd(fPsl, fPsl->blockCount-1);
        } else {
            return pslQEnd(fPsl, fIBlk);
        }
    }
    int getTBlockEnd() const {
        if (atEnd()) {
            return pslTEnd(fPsl, fPsl->blockCount-1);
        } else {
            return pslTEnd(fPsl, fIBlk);
        }
    }

    /* get offset into current block */
    int getBlockOff() const {
        return fOff;
    }
    
    /* space left in current block */
    int getBlockLeft() const {
        return getTBlockEnd() - getTPos();
    }

    /* advance by the specified amount, returning a new cursor.  If it moved
     * onto the next block, it must move to the exact beginning. If reached
     * the end, will return a cursor where atEnd() is TRUE */
    PslCursor advance(unsigned amount) const {
        assertInRange();
        assert(amount <= getBlockLeft());
        if (amount < getBlockLeft()) {
            return PslCursor(fPsl, fIBlk, fOff+amount);  // same block
        } else if (fIBlk < fPsl->blockCount-1) {
            return PslCursor(fPsl, fIBlk+1, 0); // next block
        } else {
            return PslCursor(fPsl, fPsl->blockCount, 0); // EOF
        }
    }

    /* convert to a string for debugging purposes */
    string toString() const {
        return ::toString(getQPos()) + ".." + ::toString(getQBlockEnd()) + " <> "
            + ::toString(getTPos()) + ".." + ::toString(getTBlockEnd())
            + " [" + ::toString(getBlockLeft()) + "]";
    }
};


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
 * Map one part of an exon feature.  Cursors are updated
 */
void GeneMapper::mapExonPart(const GxfFeature* exon,
                             PslCursor& srcPslCursor,
                             PslCursor& mappedPslCursor,
                             ostream& outFh) const {
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
        outFh << "    overlap: " << amount << endl;
        srcPslCursor = srcPslCursor.advance(amount);
        mappedPslCursor = mappedPslCursor.advance(amount);
    }
}

/*
 * Map an exon feature.  Cursors are updated
 */
void GeneMapper::mapExon(const GxfFeature* exon,
                         PslCursor& srcPslCursor,
                         PslCursor& mappedPslCursor,
                         ostream& outFh) const {
    assert((exon->fEnd-exon->fStart)+1 == srcPslCursor.getBlockLeft());
    outFh << "mappingExon: "  << srcPslCursor.toString() << " = " << mappedPslCursor.toString() << "\t" << exon->toString() << endl;

    // note that source blocks can be merged in mapped block, but we don't merge
    // features
    int srcPslExonQEnd = srcPslCursor.getQBlockEnd();
    while ((srcPslCursor.getQPos() < srcPslExonQEnd) && (not mappedPslCursor.atEnd())) {
        mapExonPart(exon, srcPslCursor, mappedPslCursor, outFh);
    }
    // FIXME: handled un mapped at end
    if (srcPslCursor.getQPos() < srcPslExonQEnd) {
        outFh << "   exonLeftover: " << srcPslCursor.toString() << " = " << mappedPslCursor.toString() << endl;
    }
}

/*
 * Map exons of a transcript that has at least something mapped.
 */
void GeneMapper::mapExons(const GxfFeatureVector& exons,
                          PslMapping *pslMapping,
                          ostream& outFh) const {
    PslCursor srcPslCursor(pslMapping->fSrcPsl);
    PslCursor mappedPslCursor(pslMapping->fMappedPsls[0]);
    for (int iExon = 0; iExon < exons.size(); iExon++) {
        mapExon(exons[iExon], srcPslCursor, mappedPslCursor, outFh);
    }
}

/*
 * map and output one transcripts annotations.
 */
void GeneMapper::processTranscript(const GxfFeatureNode* transcriptTree,
                                   ostream& outFh) const {
    const string& qName = transcriptTree->fFeature->getAttr("transcript_id")->fVal;
    outFh << "====" << qName  << "====" <<endl;
    GxfFeatureVector exons = getExons(transcriptTree);
    PslMapping *pslMapping =  fFeatureTransMap->mapFeatures(qName, exons);
    if (pslMapping->fMappedPsls.size() > 0) {
        mapExons(exons, pslMapping, outFh);
    } else {
        outFh << "not mapped: " << qName << endl;
    }
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
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
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
            // FIXME: outFh << gxfRecord->toString() << endl;
            delete gxfRecord;
        }
    }
}
