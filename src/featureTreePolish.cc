/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"

/* renumber ab exon */
void FeatureTreePolish::renumberExon(GxfFeature* exon,
                                     int exonNum) {
    exon->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
void FeatureTreePolish::renumberExons(FeatureNode* transcriptTree) {
    assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
    int exonNum = 1;
    // exons are always in genomic order.
    for (int i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            renumberExon(transcriptTree->fChildren[i]->fFeature, exonNum++);
        }
    }
}

/* renumber all exons in a gene */
void FeatureTreePolish::renumberGeneExons(FeatureNode* geneTreeRoot) {
    for (int i = 0; i < geneTreeRoot->fChildren.size(); i++) {
        renumberExons(geneTreeRoot->fChildren[i]);
    }
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* geneTreeRoot) {
    renumberGeneExons(geneTreeRoot);
}


