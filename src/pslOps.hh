/*
 * PSL operations and support classes.
 */
#ifndef pslOps_hh
#define pslOps_hh
#include "jkinclude.hh"
#include "typeOps.hh"
#include <vector>

/*
 * Vector of psls; doesn't own PSLs.
 */
class PslVector: public vector<struct psl*> {
    public:
    /* free all PSLs in the vector */
    void free() {
        for (int i = 0; i < size(); i++) {
            pslFree(&((*this)[i]));
        }
        clear();
    }
};

/*
 * convert a PSL to a string for debuging purposes.
 */
string pslToString(struct psl* psl);

/*
 * Convert a PSL block to a string
 */
string pslBlockToString(struct psl* psl, int blkIdx);

/* return strand as stored in psl, converting implicit `\0' to `+' */
static inline char normStrand(char strand) {
    return (strand == '\0') ? '+' : strand;
}

#if 0 // not used, may not work right
/* return query start for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslQStartStrand(struct psl *psl, int blkIdx, char strand) {
    if (psl->strand[0] == normStrand(strand)) {
        return psl->qStarts[blkIdx];
    } else {
        return psl->qSize - pslQEnd(psl, blkIdx);
    }
}

/* return query end for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslQEndStrand(struct psl *psl, int blkIdx, char strand) {
    if (psl->strand[0] == normStrand(strand)) {
        return pslQEnd(psl, blkIdx);
    } else {
        return psl->qSize - pslQStart(psl, blkIdx);
    }
}

/* return target start for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslTStartStrand(struct psl *psl, int blkIdx, char strand) {
    if (normStrand(psl->strand[1]) == normStrand(strand)) {
        return psl->tStarts[blkIdx];
    } else {
        return psl->tSize - pslTEnd(psl, blkIdx);
    }
}

/* return target end for the given block, mapped to specified strand,
 * which can be `\0' for `+' */
static inline unsigned pslTEndStrand(struct psl *psl, int blkIdx, char strand) {
    if (normStrand(psl->strand[1]) == normStrand(strand)) {
        return pslTEnd(psl, blkIdx);
    } else {
        return psl->tSize - pslTStart(psl, blkIdx);
    }
}
#endif

/* is the query fully mapped */
static inline bool pslQueryFullyMapped(struct psl* psl) {
    // look for any query gaps
    if (psl->qStart > 0) {
        return false;
    }
    for (int iBlk = 1; iBlk < psl->blockCount; iBlk++) {
        if (pslQStart(psl, iBlk) != pslQEnd(psl, iBlk-1)) {
            return false;
        }
    }
    if (psl->qEnd < psl->qSize) {
        return false;
    }
    return true;
}

/* add a psl block to a psl being constructed, updating counts */
static inline void pslAddBlock(struct psl* psl, unsigned qStart, unsigned tStart, int blockSize) {
    int iBlk = psl->blockCount;
    psl->qStarts[iBlk] = qStart;
    psl->tStarts[iBlk] = tStart;
    psl->blockSizes[iBlk] = blockSize;
    psl->match += blockSize;
    if (iBlk > 0) {
        if (pslQStart(psl, iBlk) > pslQEnd(psl, iBlk-1)) {
            psl->qNumInsert++;
            psl->qBaseInsert += pslQStart(psl, iBlk)-pslQEnd(psl, iBlk-1);
        }
        if (pslTStart(psl, iBlk) > pslTEnd(psl, iBlk-1)) {
            psl->tNumInsert++;
            psl->tBaseInsert += pslTStart(psl, iBlk)-pslTEnd(psl, iBlk-1);
        }
    }
    psl->blockCount++;
}

#endif
