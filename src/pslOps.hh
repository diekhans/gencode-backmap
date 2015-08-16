/*
 * PSL operations and support classes.
 */
#ifndef pslOps_hh
#define pslOps_hh
#include "jkinclude.hh"
#include "typeOps.hh"


/*
 * convert a PSL to a string for debuging purposes.
 */
string pslToString(struct psl* psl);

/*
 * Convert a PSL block to a string
 */
string pslBlockToString(struct psl* psl, int blkIdx);

/* return query start for the given block */
static inline unsigned pslQStart(struct psl *psl, int blkIdx) {
    return psl->qStarts[blkIdx];
}

/* return target start for the given block */
static inline unsigned pslTStart(struct psl *psl, int blkIdx) {
    return psl->tStarts[blkIdx];
}

/* return strand as stored in psl, converting implicit `\0' to `+' */
static inline char normStrand(char strand) {
    return (strand == '\0') ? '+' : strand;
}

/* return query start for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslQStartStrand(struct psl *psl, int blkIdx, char strand) {
    if (psl->strand[0] == normStrand(strand)) {
        return psl->qStarts[blkIdx];
    } else {
        return psl->qSize - psl->qStarts[blkIdx];
    }
}

/* return query end for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslQEndStrand(struct psl *psl, int blkIdx, char strand) {
    if (psl->strand[0] == normStrand(strand)) {
        return pslQEnd(psl, blkIdx);
    } else {
        return psl->qSize - pslQEnd(psl, blkIdx);
    }
}

/* return target start for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslTStartStrand(struct psl *psl, int blkIdx, char strand) {
    if (normStrand(psl->strand[1]) == normStrand(strand)) {
        return psl->tStarts[blkIdx];
    } else {
        return psl->tSize - psl->tStarts[blkIdx];
    }
}

/* return target end for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslTEndStrand(struct psl *psl, int blkIdx, char strand) {
    if (normStrand(psl->strand[1]) == normStrand(strand)) {
        return pslTEnd(psl, blkIdx);
    } else {
        return psl->tSize - pslTEnd(psl, blkIdx);
    }
}

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
    PslCursor(const PslCursor& src):
        fPsl(src.fPsl), fIBlk(src.fIBlk), fOff(src.fOff) {
    }
    struct psl* getPsl() const {
        return fPsl;
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

    /* accessors of current position, strand adjusted */
    int getQPosStrand(char strand) const {
        if (pslQStrand(fPsl) == strand) {
            return getQPos();
        } else {
            return fPsl->qSize - getQBlockEnd();
        }
    }
    int getTPosStrand(char strand) const {
        if (pslTStrand(fPsl) == strand) {
            return getTPos();
        } else {
            return fPsl->tSize - getTBlockEnd();
        }
    }

    /* accessors for current block end, strand adjusted */
    int getQBlockEndStrand(char strand) const {
        if (strand == pslQStrand(fPsl)) {
            return getQBlockEnd();
        } else {
            return fPsl->qSize - getQPos();
        }
    }
    int getTBlockEndStrand(char strand) const {
        if (strand == pslTStrand(fPsl)) {
            return getTBlockEnd();
        } else {
            return fPsl->tSize - getTPos();
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

#endif
