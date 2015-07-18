/*
 * transmap projection of annotations
 */
#include "transMapper.hh"
#include "jkcommon.hh"
extern "C" {
#define hash jkhash
#define new jknew
#include "pslTransMap.h"
#include "psl.h"
#include "genomeRangeTree.h"
#include "chain.h"
#include "dnautil.h"
#undef hash
#undef new
}

/* slCat that reverses parameter order, as the first list in rangeTreeAddVal
 * mergeVals function tends to be larger in degenerate cases of a huge number
 * of chains */
static void *slCatReversed(void *va, void *vb) {
    return slCat(vb, va);
}

/* add a map align object to the genomeRangeTree */
void TransMapper::mapAlnsAdd(struct psl *mapPsl) {
    genomeRangeTreeAddVal(fMapAlns, mapPsl->qName, mapPsl->qStart, mapPsl->qEnd, mapPsl, slCatReversed);
}

/* convert a chain to a psl, ignoring match counts, etc */
struct psl* TransMapper::chainToPsl(struct chain *ch,
                                    bool swapMap) {
    struct psl *psl;
    struct cBlock *cBlk;
    int iBlk;
    int qStart = ch->qStart, qEnd = ch->qEnd;
    char strand[2] = {ch->qStrand, '\0'};
    if (ch->qStrand == '-') {
        reverseIntRange(&qStart, &qEnd, ch->qSize);
    }
    psl = pslNew(ch->qName, ch->qSize, qStart, qEnd,
                 ch->tName, ch->tSize, ch->tStart, ch->tEnd,
                 strand, slCount(ch->blockList), 0);
    for (cBlk = ch->blockList, iBlk = 0; cBlk != NULL; cBlk = cBlk->next, iBlk++) {
        psl->blockSizes[iBlk] = (cBlk->tEnd - cBlk->tStart);
        psl->qStarts[iBlk] = cBlk->qStart;
        psl->tStarts[iBlk] = cBlk->tStart;
        psl->match += psl->blockSizes[iBlk];
    }
    psl->blockCount = iBlk;
    if (swapMap)
        pslSwap(psl, FALSE);
    return psl;
}

/* read a chain file, convert to mapAln object and genomeRangeTree by query locations. */
void TransMapper::loadMapChains(const string& chainFile,
                                bool swapMap) {
    struct genomeRangeTree* mapAlns = genomeRangeTreeNew();
    struct chain *ch;
    struct lineFile *chLf = lineFileOpen(const_cast<char*>(chainFile.c_str()), TRUE);
    while ((ch = chainRead(chLf)) != NULL) {
        mapAlnsAdd(chainToPsl(ch, swapMap));
        chainFree(&ch);
    }
    lineFileClose(&chLf);
}

/* constructor, loading chains */
TransMapper::TransMapper(const string& chainFile,
                         bool swapMap):
    fMapAlns(genomeRangeTreeNew()) {
    loadMapChains(chainFile, swapMap);
}

/* map one pair of query and target PSL */
struct psl* TransMapper::mapPslPair(struct psl *inPsl, struct psl *mapPsl) {
    if (inPsl->tSize != mapPsl->qSize)
        errAbort(const_cast<char*>("Error: inPsl %s tSize (%d) != mapping alignment %s qSize (%d) (perhaps you need to specify -swapMap?)"),
                 inPsl->tName, inPsl->tSize, mapPsl->qName, mapPsl->qSize);
    return pslTransMap(pslTransMapNoOpts, inPsl, mapPsl);
}

/* Map a single input PSL and return a list of resulting mappings */
struct psl* TransMapper::mapPsl(struct psl* inPsl) {
    struct psl* mappedPsls = NULL;
    struct range *overMapAlnNodes = genomeRangeTreeAllOverlapping(fMapAlns, inPsl->tName, inPsl->tStart, inPsl->tEnd);
    for (struct range *overMapAlnNode = overMapAlnNodes; overMapAlnNode != NULL; overMapAlnNode = overMapAlnNode->next) {
        for (struct psl *overMapPsl = static_cast<struct psl*>(overMapAlnNode->val); overMapPsl != NULL; overMapPsl = overMapPsl->next) {
            struct psl* mappedPsl = mapPslPair(inPsl, overMapPsl);
            if (mappedPsl != NULL) {
                slAddHead(&mappedPsls, mappedPsl);
            }
        }
    }
    return mappedPsls;
}
