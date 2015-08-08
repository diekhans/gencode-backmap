#include "typeOps.hh"
#include <stdexcept>
#include <errno.h>
#include <stdlib.h>
#include "jkinclude.hh"

const string whitespace = " \t\n\r\f\v";


/*
 * Convert a string to an int.
 */
int stringToInt(const string& str,
                bool* isOk,
                int base) {
    const char* cstr = str.c_str();
    char *endPtr;
    errno = 0;
    long lnum = strtol(cstr, &endPtr, base);
    if ((endPtr == cstr) || (*endPtr != '\0')) {
        if (isOk != NULL) {
            *isOk = false;
            return 0;
        } else {
            throw invalid_argument("Invalid integer \"" + string(str) + "\"");
        }
    }
     
    int num = (int)lnum;
    if ((errno != 0) || ((long)num != lnum)) {
        if (isOk != NULL) {
            *isOk = false;
            return 0;
        } else {
            throw invalid_argument("Integer out of range \"" + string(str) + "\"");
        }
    }
    if (isOk != NULL) {
        *isOk = true;
    }
    return num;
}

/*
 * Split a string into a vector of string given a separator character.
 */
StringVector stringSplit(const string& str,
                         char separator) {
    StringVector strs;
    
    int prevIdx = 0;
    int sepIdx;
    while ((sepIdx = (int)str.find_first_of(separator, prevIdx)) >= 0) {
        strs.push_back(str.substr(prevIdx, sepIdx-prevIdx));
        prevIdx = sepIdx+1;
    }
    strs.push_back(str.substr(prevIdx));

    return strs;
}

/*
 * Convert an integer to a string.
 */
string toString(int num) {
    char buf[64];
    sprintf(buf, "%d", num);
    return string(buf);
}

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
