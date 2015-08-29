/*
 * transmap projection of annotations
 */
#include "transMap.hh"
#include "jkinclude.hh"
#include "typeOps.hh"

/* slCat that reverses parameter order, as the first list in rangeTreeAddVal
 * mergeVals function tends to be larger in degenerate cases of a huge number
 * of chains */
static void *slCatReversed(void *va, void *vb) {
    return slCat(vb, va);
}

/* add a query or target size if it doesn't already exist */
void TransMap::addSeqSize(const string& seqName,
                          int seqSize,
                          SizeMap& sizeMap) {
    if (sizeMap.find(seqName) == sizeMap.end()) {
        sizeMap[seqName] = seqSize;
    }
}

/* add a map align object to the genomeRangeTree */
void TransMap::mapAlnsAdd(struct psl *mapPsl) {
    genomeRangeTreeAddVal(fMapAlns, mapPsl->qName, mapPsl->qStart, mapPsl->qEnd, mapPsl, slCatReversed);
    addSeqSize(mapPsl->qName, mapPsl->qSize, fQuerySizes);
    addSeqSize(mapPsl->tName, mapPsl->tSize, fTargetSizes);
}

/* convert a chain to a psl, ignoring match counts, etc */
struct psl* TransMap::chainToPsl(struct chain *ch,
                                    bool swapMap) {
    int qStart = ch->qStart, qEnd = ch->qEnd;
    char strand[2] = {ch->qStrand, '\0'};
    if (ch->qStrand == '-') {
        reverseIntRange(&qStart, &qEnd, ch->qSize);
    }
    struct psl* psl = pslNew(ch->qName, ch->qSize, qStart, qEnd,
                             ch->tName, ch->tSize, ch->tStart, ch->tEnd,
                             strand, slCount(ch->blockList), 0);
    int iBlk = 0;
    for (struct cBlock *cBlk = ch->blockList; cBlk != NULL; cBlk = cBlk->next, iBlk++) {
        pslAddBlock(psl, cBlk->qStart, cBlk->tStart, (cBlk->tEnd - cBlk->tStart));
    }
    if (swapMap)
        pslSwap(psl, FALSE);
    return psl;
}

/* read a chain file, convert to mapAln object and genomeRangeTree by query locations. */
void TransMap::loadMapChains(const string& chainFile,
                             bool swapMap) {
    struct chain *ch;
    struct lineFile *chLf = lineFileOpen(toCharStr(chainFile), TRUE);
    while ((ch = chainRead(chLf)) != NULL) {
        mapAlnsAdd(chainToPsl(ch, swapMap));
        chainFree(&ch);
    }
    lineFileClose(&chLf);
}

/* constructor */
TransMap::TransMap():
    fMapAlns(genomeRangeTreeNew()) {
}

/* destructor */
TransMap::~TransMap() {
    struct hashCookie cookie = hashFirst(fMapAlns->jkhash);
    struct hashEl* chromEl;
    while ((chromEl = hashNext(&cookie)) != NULL) {
        for (struct range* range = genomeRangeTreeList(fMapAlns, chromEl->name);
             range != NULL; range = range->next) {
            struct psl* head = static_cast<struct psl*>(range->val);
            pslFreeList(&head);
        }
    }
    
    genomeRangeTreeFree(&fMapAlns);
}



/* map one pair of query and target PSL */
struct psl* TransMap::mapPslPair(struct psl *inPsl, struct psl *mapPsl) const {
    if (inPsl->tSize != mapPsl->qSize)
        errAbort(toCharStr("Error: inPsl %s tSize (%d) != mapping alignment %s qSize (%d) (perhaps you need to specify -swapMap?)"),
                 inPsl->tName, inPsl->tSize, mapPsl->qName, mapPsl->qSize);
    return pslTransMap(pslTransMapKeepTrans, inPsl, mapPsl);
}

/* Map a single input PSL and return a list of resulting mappings.
 * Keep PSL in the same order, even if it creates a `-' on the target. */
PslVector TransMap::mapPsl(struct psl* inPsl) const {
    PslVector mappedPsls;
    struct range *overMapAlnNodes = genomeRangeTreeAllOverlapping(fMapAlns, inPsl->tName, inPsl->tStart, inPsl->tEnd);
    for (struct range *overMapAlnNode = overMapAlnNodes; overMapAlnNode != NULL; overMapAlnNode = overMapAlnNode->next) {
        for (struct psl *overMapPsl = static_cast<struct psl*>(overMapAlnNode->val); overMapPsl != NULL; overMapPsl = overMapPsl->next) {
            struct psl* mappedPsl = mapPslPair(inPsl, overMapPsl);
            if (mappedPsl != NULL) {
                mappedPsls.push_back(mappedPsl);
            }
        }
    }
    return mappedPsls;
}

/* factory from a chain file */
TransMap* TransMap::factoryFromChainFile(const string& chainFile,
                                         bool swapMap) {
    TransMap* transMap = new TransMap();
    transMap->loadMapChains(chainFile, swapMap);
    return transMap;
}

/* factory from a list of psls */
TransMap* TransMap::factoryFromPsls(struct psl* psls,
                                    bool swapMap) {
    TransMap* transMap = new TransMap();
    for (struct psl* psl = psls; psl != NULL; psl = psl->next) {
        struct psl* pslCp = pslClone(psl);
        if (swapMap)
            pslSwap(pslCp, FALSE);
        transMap->mapAlnsAdd(pslCp);
    }
    return transMap;
}
    
/* factory from a psl file */
TransMap* TransMap::factoryFromPslFile(const string& pslFile,
                                       bool swapMap) {
    struct psl* psls = pslLoadAll(toCharStr(pslFile));
    TransMap* transMap = factoryFromPsls(psls, swapMap);
    pslFreeList(&psls);
    return transMap;
}
