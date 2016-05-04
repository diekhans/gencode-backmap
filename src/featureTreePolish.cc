/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>

/* renumber ab exon */
void FeatureTreePolish::renumberExon(FeatureNode* exonNode,
                                     int exonNum,
                                     ExonIdExonMap& exonIdExonMap) {
    int oldExonNum = stringToInt(exonNode->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR)->getVal());
    exonIdExonMap[oldExonNum].push_back(exonNode);
    exonNode->fFeature->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
void FeatureTreePolish::renumberExons(FeatureNode* transcriptTree,
                                      ExonIdExonMap& exonIdExonMap) {
    int exonNum = 1;
    // exons are always in genomic order.
    for (int i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            renumberExon(transcriptTree->fChildren[i], exonNum++, exonIdExonMap);
        }
    }
}

/* find the new exon containing a feature given the old exon number */
FeatureNode* FeatureTreePolish::findNewExon(FeatureNode* featureNode,
                                            int oldExonNum,
                                            ExonIdExonMap& exonIdExonMap) {
    for (int i = 0; i < exonIdExonMap[oldExonNum].size(); i++) {
        if (featureNode->fFeature->overlaps(exonIdExonMap[oldExonNum][i]->fFeature)) {
            return exonIdExonMap[oldExonNum][i];
        }
    }
    throw logic_error("renumberOtherFeature: lost exon");
}

/* change a non-exon feature to their exon numbers */
void FeatureTreePolish::renumberOtherFeature(FeatureNode* featureNode,
                                             ExonIdExonMap& exonIdExonMap) {
    const AttrVal* exonNumAttr = featureNode->fFeature->getAttrs().find(GxfFeature::EXON_NUMBER_ATTR);;
    if (exonNumAttr != NULL) {
        FeatureNode *newExon = findNewExon(featureNode, stringToInt(exonNumAttr->getVal()), exonIdExonMap);
        const AttrVal* newExonNumAttr = newExon->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR);
        featureNode->fFeature->getAttrs().update(*newExonNumAttr);
    }
}

/* recursively change non-exons features to match the news exon numbers */
void FeatureTreePolish::renumberOtherFeatures(FeatureNode* featureNode,
                                              ExonIdExonMap& exonIdExonMap) {
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        if (featureNode->fChildren[i]->fFeature->fType != GxfFeature::EXON) {
            renumberOtherFeature(featureNode->fChildren[i], exonIdExonMap);
        }
        renumberOtherFeatures(featureNode->fChildren[i], exonIdExonMap);
    }
}

/* renumber all features in a transcript */
void FeatureTreePolish::renumberTranscript(FeatureNode* transcriptTree) {
    assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
    ExonIdExonMap exonIdExonMap;
    renumberExons(transcriptTree, exonIdExonMap);
    renumberOtherFeatures(transcriptTree, exonIdExonMap);
}

/* renumber all exons in a gene */
void FeatureTreePolish::renumberGeneExons(FeatureNode* geneTreeRoot) {
    for (int i = 0; i < geneTreeRoot->fChildren.size(); i++) {
        renumberTranscript(geneTreeRoot->fChildren[i]);
    }
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* geneTreeRoot) {
    renumberGeneExons(geneTreeRoot);
}


