
#include "geneMapper.hh"
#include "jkcommon.hh"
extern "C" {
#define hash jkhash
#define new jknew
#include "psl.h"
#include "dnautil.h"
#undef hash
#undef new
}
#include "transMapper.hh"
#include "gxf.hh"
#include <algorithm>
#include <fstream>

/* constructor, sort mapped PSLs */
PslMapping::PslMapping(struct psl* srcPsl,
                       PslVector& mappedPsls):
    fSrcPsl(srcPsl),
    fMappedPsls(mappedPsls),
    fScore(0) {
    sortMappedPsls();
    fScore = calcPslMappingScore(srcPsl, mappedPsls[0]); // best score
}

/* free up psls */
PslMapping::~PslMapping() {
    pslFree(&fSrcPsl);
    for (size_t i = 0; i < fMappedPsls.size(); i++) {
        pslFree(&(fMappedPsls[i]));
    }
}

/* calculate the number of aligned bases */
int PslMapping::numPslAligned(struct psl* psl) {
    return (psl->match + psl->repMatch + psl->misMatch);
}

/* Compute a mapping score between a src and mapped psl.  A perfect mapping is
 * a zero score.  Extra inserts count against the score.
 */
int PslMapping::calcPslMappingScore(struct psl* srcPsl,
                                    struct psl* mappedPsl) {
    return (numPslAligned(srcPsl) - numPslAligned(mappedPsl))
        + abs(int(srcPsl->qNumInsert)-int(mappedPsl->qNumInsert))
        + abs(int(srcPsl->tNumInsert)-int(mappedPsl->tNumInsert));
}

/* comparison functor based on lowest store. */
class ScoreCmp {
    public:
    struct psl* fSrcPsl;
    
    ScoreCmp(struct psl* srcPsl):
        fSrcPsl(srcPsl) {
    }

    int operator()(struct psl* mappedPsl1,
                   struct psl* mappedPsl2) const {
        return -(PslMapping::calcPslMappingScore(fSrcPsl, mappedPsl1)
                 - PslMapping::calcPslMappingScore(fSrcPsl, mappedPsl2));
    }
};
/* sort with best (lowest score) first */
void PslMapping::sortMappedPsls() {
    ScoreCmp scoreCmp(fSrcPsl);
    sort(fMappedPsls.begin(), fMappedPsls.end(), scoreCmp);
}


/** compute total size of features */
int GeneMapper::sumFeatureSizes(const GxfFeatureVector& features) {
    int sz = 0;
    for (size_t i = 0; i < features.size(); i++) {
        sz += (features[i]->fEnd - features[i]->fStart)+1;
    }
    return sz;
}
/* create a psl from a list of features. assumes features are sorter in
 * ascending order */
struct psl* GeneMapper::featuresToPsl(const string& qName,
                                      const GxfFeatureVector& exons) {
    // this does a [1..n] to [0..n) conversion
    int qSize = sumFeatureSizes(exons);
    int tSize = fTransMapper->getQuerySeqSize(exons[0]->fSeqid); // target is mapping query
    struct psl* psl = pslNew(toCharStr(qName), qSize, 0, qSize,
                             toCharStr(exons[0]->fSeqid), tSize, exons[0]->fStart-1, exons[exons.size()-1]->fEnd,
                             toCharStr(exons[0]->fStrand), exons.size(), 0);
    int qStart = 0;
    for (size_t iBlk = 0; iBlk < exons.size(); iBlk++) {
        const GxfFeature* feature = exons[iBlk];
        psl->blockSizes[iBlk] = (feature->fEnd - feature->fStart)+1;
        psl->qStarts[iBlk] = qStart;
        psl->tStarts[iBlk] = feature->fStart-1;
        psl->match += psl->blockSizes[iBlk];
        qStart += psl->blockSizes[iBlk];
        if ((iBlk > 0) && (psl->tStarts[iBlk] > (psl->tStarts[iBlk-1]+psl->blockSizes[iBlk]))) {
            throw invalid_argument("out of order blocks: " + feature->toString());
        }
    }
    return psl;
}

/* get exon features */
GxfFeatureVector GeneMapper::getExons(const GxfFeatureNode* transcriptTree) {
    GxfFeatureVector exons;
    for (size_t i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            exons.push_back(transcriptTree->fChildren[i]->fFeature);
        }
    }
    return exons;
}

/* create a psl for exons of a transcript */
struct psl* GeneMapper::transcriptExonsToPsl(const GxfFeatureNode* transcriptTree) {
    const string& qName = transcriptTree->fFeature->getAttr("transcript_id")->fVal;
    const GxfFeatureVector exons = getExons(transcriptTree);
    return featuresToPsl(qName, exons);
}

/* 
 * Map a transcript's exons, returning object that has mappings and scoring
 */
PslMapping* GeneMapper::mapTranscriptExons(const GxfFeatureNode* transcriptTree) {
    struct psl* transcriptPsl = transcriptExonsToPsl(transcriptTree);
    PslVector mappedTranscriptPsls = fTransMapper->mapPsl(transcriptPsl);
    return new PslMapping(transcriptPsl, mappedTranscriptPsls);
}


/* Return a list records, moving from gxfRecords vector  */
void GeneMapper::queueRecords(GxfParser *gxfParser,
                              GxfRecordVector& gxfRecords) {
    for (size_t i = 0; i < gxfRecords.size(); i++) {
        gxfParser->push(gxfRecords[i]);
    }
    gxfRecords.clear();
}

/*
 * Find the parent for GFF3.
 */
GxfFeatureNode* GeneMapper::findGff3Parent(GxfFeatureNode* geneTreeLeaf,
                                           const GxfFeature* gxfFeature) {
    const string& parentId = gxfFeature->getAttrValue("Parent");
    GxfFeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->getAttrValue("ID") != parentId)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node not found: " + parentId);
    }
    return parent;
}

/*
 * Process a FF3 record for a gene, which uses the explicit tree.
 * Return the new leaf node.
 */
GxfFeatureNode* GeneMapper::loadGff3GeneRecord(const GxfFeature* gxfFeature,
                                               GxfFeatureNode* geneTreeLeaf) {
    GxfFeatureNode* parent = findGff3Parent(geneTreeLeaf, gxfFeature);
    parent->addChild(new GxfFeatureNode(gxfFeature));
    return parent;
}

/* Get the desired type the parent of GTF feature
 * WARNING: this assumed that the hierarchy is
 * gene->transcript->{everything else}
 * FIXME: this is what GFF3 does, which might not be right.
 */
const string& GeneMapper::getGtfParentType(const string& featureType) {
    assert(featureType != GxfFeature::GENE);
    if (featureType == GxfFeature::TRANSCRIPT) {
        return GxfFeature::GENE;
    } else {
        return GxfFeature::TRANSCRIPT;
    }
}

/*
 * Find the parent for a GTF record.  This is painful guess based on the
 * GENCODE file order and know how GENCODE is structures.
 */
GxfFeatureNode* GeneMapper::findGtfParent(GxfFeatureNode* geneTreeLeaf,
                                          const GxfFeature* gxfFeature) {
    const string& parentType = getGtfParentType(gxfFeature->fType);
    GxfFeatureNode* parent = geneTreeLeaf;
    while ((parent != NULL) and (parent->fFeature->fType != parentType)) {
        parent = parent->fParent;
    }
    if (parent == NULL) {
        throw invalid_argument("parent node not found: " + parentType);
    }
    return parent;
}

/*
 * Process a GTF record for a gene, which uses knowledge of
 * the GENCODE structure to reproduce the hierarchy.
 * Return the new leaf node.
 */
GxfFeatureNode* GeneMapper::loadGtfGeneRecord(const GxfFeature* gxfFeature,
                                              GxfFeatureNode* geneTreeLeaf) {
    GxfFeatureNode* parent = findGtfParent(geneTreeLeaf, gxfFeature);
    parent->addChild(new GxfFeatureNode(gxfFeature));
    return parent;
}

/*
 * Process a GxfRecord for a gene, return False if no more for this gene.
 */
bool GeneMapper::loadGeneRecord(GxfParser *gxfParser,
                                const GxfRecord* gxfRecord,
                                GxfFeatureNode* geneTreeRoot,
                                GxfFeatureNode*& geneTreeLeaf,
                                GxfRecordVector& queuedRecords) {
    if (instanceOf(gxfRecord, GxfLine)) {
        queuedRecords.push_back(gxfRecord);
        return true;
    } else {
        const GxfFeature* gxfFeature = dynamic_cast<const GxfFeature*>(gxfRecord);
        if (gxfFeature->fType == GxfFeature::GENE) {
            queuedRecords.push_back(gxfRecord); // next gene
            return false;
        } else {
            if (gxfParser->getGxfFormat() == GFF3_FORMAT) {
                geneTreeLeaf = loadGff3GeneRecord(gxfFeature, geneTreeLeaf);
            } else {
                geneTreeLeaf = loadGtfGeneRecord(gxfFeature, geneTreeLeaf);
            }
            return true;
        }
    }
}

/*
 * When a new gene record has been reached, load it.  Return non-feature and
 * the next gene to the queue to process.  This will causes comments in the
 * middle of genes to be moved to the end, but GENCODE doesn't do this.  This
 * whole thing is annoying due to the lack of explicit structure in GTF.
 */
GxfFeatureNode* GeneMapper::loadGene(GxfParser *gxfParser,
                                     const GxfFeature* geneFeature) {
    assert(geneFeature->fType == GxfFeature::GENE);

    GxfFeatureNode* geneTreeRoot = new GxfFeatureNode(geneFeature);
    GxfFeatureNode* geneTreeLeaf = geneTreeRoot;  // were we are currently working
    GxfRecordVector queuedRecords;
    const GxfRecord* gxfRecord = NULL;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (not loadGeneRecord(gxfParser, gxfRecord, geneTreeRoot, geneTreeLeaf, queuedRecords)) {
            break;
        }
    }
    queueRecords(gxfParser, queuedRecords);
    return geneTreeRoot;
}

/* Map a GFF3/GTF */
void GeneMapper::mapGxf(GxfParser *gxfParser,
                        ostream& outFh) {
    const GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        if (instanceOf(gxfRecord, GxfFeature)) {
            const GxfFeature* geneFeature = dynamic_cast<const GxfFeature*>(gxfRecord);
            GxfFeatureNode* geneTree = loadGene(gxfParser, geneFeature);
            delete geneTree;
        } else {
            // outGxf << gxfRecord->toString() << endl;
            delete gxfRecord;
        }
    }
}
