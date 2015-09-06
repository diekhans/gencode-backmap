/*
 * Cursor used to traverse PSLs.
 */
#ifndef pslCursor_hh
#define pslCursor_hh
#include "jkinclude.hh"
#include "typeOps.hh"

/*
 * Cursor into a PSL.  Tracks position in an alignment.
 */
class PslCursor {
    private:
    struct psl* fPsl;
    int fBlkIdx;   // index of current block, set to blockCount if reached end
    int fBlkOff;   // offset in current block;

    /* check if in range */
    inline void assertInRange() const {
        assert(fBlkIdx < fPsl->blockCount);
        assert(fBlkOff < fPsl->blockSizes[fBlkIdx]);
    }
    
    public:
    PslCursor(struct psl* psl, int iBlk=0, int off=0):
        fPsl(psl), fBlkIdx(iBlk), fBlkOff(off) {
    }
    PslCursor(const PslCursor& src):
        fPsl(src.fPsl), fBlkIdx(src.fBlkIdx), fBlkOff(src.fBlkOff) {
    }
    struct psl* getPsl() const {
        return fPsl;
    }
    
    /* have we reached the end of the psl */
    bool atEnd() const {
        return fBlkIdx >= fPsl->blockCount;
    }

    /* accessors of current position */
    int getQPos() const {
        if (atEnd()) {
            return pslQEnd(fPsl, fPsl->blockCount-1); 
        } else {
            assertInRange();
            return pslQStart(fPsl, fBlkIdx)+fBlkOff;
         }
    }
    int getTPos() const {
        if (atEnd()) {
            return pslTEnd(fPsl, fPsl->blockCount-1); 
        } else {
            assertInRange();
            return pslTStart(fPsl, fBlkIdx)+fBlkOff;
        }
    }

    /* accessors of start position */
    int getQStart() const {
        return pslQStart(fPsl, 0);
    }

    /* accessors for current block end */
    int getQBlockEnd() const {
        if (atEnd()) {
            return pslQEnd(fPsl, fPsl->blockCount-1);
        } else {
            return pslQEnd(fPsl, fBlkIdx);
        }
    }
    int getTBlockEnd() const {
        if (atEnd()) {
            return pslTEnd(fPsl, fPsl->blockCount-1);
        } else {
            return pslTEnd(fPsl, fBlkIdx);
        }
    }

    /* accessor for range at current position, strand adjusted */
    void getTRangeStrand(char strand, int length, int* tStart, int* tEnd) const {
        if (pslTStrand(fPsl) == strand) {
            *tStart = getTPos();
            *tEnd = getTPos()+length;
        } else {
            *tStart = (fPsl->tSize-(getTPos()+length));
            *tEnd = (fPsl->tSize-getTPos());
        }
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
            return PslCursor(fPsl, fBlkIdx, fBlkOff+amount);  // same block
        } else if (fBlkIdx < fPsl->blockCount-1) {
            return PslCursor(fPsl, fBlkIdx+1, 0); // next block
        } else {
            return PslCursor(fPsl, fPsl->blockCount, 0); // EOF
        }
    }

    /* convert to a string for debugging purposes, if length is supplied display end to that length */
    string toString(int length=-1) const {
        int qPos = getQPos();
        int tPos = getTPos();
        int qEnd = (length < 0) ? getQBlockEnd() : (qPos+length);
        int tEnd = (length < 0) ? getTBlockEnd() : (tPos + length);
        string str = string(fPsl->qName) + ":" + charToString(pslQStrand(fPsl)) + ":"
            + ::toString(qPos) + ".." + ::toString(qEnd) + " <|> "
            + string(fPsl->tName) + ":" + charToString(pslTStrand(fPsl)) + ":"
            + ::toString(tPos) + ".." + ::toString(tEnd);
        // reverse strand too
        char tStrand = pslTStrand(fPsl) == '+' ? '-' : '+';
        char qStrand = pslQStrand(fPsl) == '+' ? '-' : '+';
        reverseIntRange(&qPos, &qEnd, fPsl->qSize);
        reverseIntRange(&tPos, &tEnd, fPsl->tSize);
        str += " (" + charToString(qStrand) + ":"
            + ::toString(qPos) + ".." + ::toString(qEnd) + " <|> "
            +  charToString(tStrand) + ":" + ::toString(tPos) + ".." + ::toString(tEnd) + ")";
        str += " [" + ::toString(fBlkIdx) + ", " + ::toString(fBlkOff) + ", " + ::toString(getBlockLeft()) + "]";
        return str;
    }
};

#endif
