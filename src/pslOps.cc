#include "pslOps.hh"


/*
 * convert an autoSql unsiged array to a commastring */
static string unsignedArrayToString(unsigned len,
                                    const unsigned* vals) {
    string valStr;
    for (unsigned i = 0; i < len; i++) {
        valStr += toString(vals[i]) + ",";
    }
    return valStr;
}

/*
 * convert a PSL to a string for debuging purposes.
 */
string pslToString(struct psl* psl) {
    return toString(psl->match) + "\t" +
        toString(psl->misMatch) + "\t" +
        toString(psl->repMatch) + "\t" +
        toString(psl->nCount) + "\t" +
        toString(psl->qNumInsert) + "\t" +
        toString(psl->qBaseInsert) + "\t" +
        toString(psl->tNumInsert) + "\t" +
        toString(psl->tBaseInsert) + "\t" +
        string(psl->strand) + "\t" +
        string(psl->qName) + "\t" +
        toString(psl->qSize) + "\t" +
        toString(psl->qStart) + "\t" +
        toString(psl->qEnd) + "\t" +
        string(psl->tName) + "\t" +
        toString(psl->tSize) + "\t" +
        toString(psl->tStart) + "\t" +
        toString(psl->tEnd) + "\t" +
        toString(psl->blockCount) + "\t" +
        unsignedArrayToString(psl->blockCount, psl->blockSizes) + "\t" +
        unsignedArrayToString(psl->blockCount, psl->qStarts) + "\t" +
        unsignedArrayToString(psl->blockCount, psl->tStarts);
}

/*
 * Convert a PSL block to a string
 */
string pslBlockToString(struct psl* psl, int blkIdx) {
    return string(psl->qName) + ":" + charToString(psl->strand[0]) + ":"
        + toString(pslQStart(psl, blkIdx)) + "-"
        + toString(pslQEnd(psl, blkIdx))
        + " => " +
        string(psl->tName) + ":" + charToString(normStrand(psl->strand[1])) + ":"
        + toString(pslTStart(psl, blkIdx)) + "-"
        + toString(pslTEnd(psl, blkIdx));
}

