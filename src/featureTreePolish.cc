/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>

/* renumber ab exon */
void FeatureTreePolish::renumberExon(FeatureNode* exonNode,
                                     int exonNum,
                                     ExonIdExonMap& exonIdExonMap) const {
    int oldExonNum = stringToInt(exonNode->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR)->getVal());
    exonIdExonMap[oldExonNum].push_back(exonNode);
    exonNode->fFeature->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
void FeatureTreePolish::renumberExons(FeatureNode* transcriptTree,
                                      ExonIdExonMap& exonIdExonMap) const {
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
                                            ExonIdExonMap& exonIdExonMap) const {
    for (int i = 0; i < exonIdExonMap[oldExonNum].size(); i++) {
        if (featureNode->fFeature->overlaps(exonIdExonMap[oldExonNum][i]->fFeature)) {
            return exonIdExonMap[oldExonNum][i];
        }
    }
    throw logic_error("renumberOtherFeature: lost exon");
}

/* change a non-exon feature to their exon numbers */
void FeatureTreePolish::renumberOtherFeature(FeatureNode* featureNode,
                                             ExonIdExonMap& exonIdExonMap) const {
    const AttrVal* exonNumAttr = featureNode->fFeature->getAttrs().find(GxfFeature::EXON_NUMBER_ATTR);;
    if (exonNumAttr != NULL) {
        FeatureNode *newExon = findNewExon(featureNode, stringToInt(exonNumAttr->getVal()), exonIdExonMap);
        const AttrVal* newExonNumAttr = newExon->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR);
        featureNode->fFeature->getAttrs().update(*newExonNumAttr);
    }
}

/* recursively change non-exons features to match the news exon numbers */
void FeatureTreePolish::renumberOtherFeatures(FeatureNode* featureNode,
                                              ExonIdExonMap& exonIdExonMap) const {
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        if (featureNode->fChildren[i]->fFeature->fType != GxfFeature::EXON) {
            renumberOtherFeature(featureNode->fChildren[i], exonIdExonMap);
        }
        renumberOtherFeatures(featureNode->fChildren[i], exonIdExonMap);
    }
}

/* renumber all features in a transcript */
void FeatureTreePolish::renumberTranscriptExons(FeatureNode* transcriptTree) const {
    assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
    ExonIdExonMap exonIdExonMap;
    renumberExons(transcriptTree, exonIdExonMap);
    renumberOtherFeatures(transcriptTree, exonIdExonMap);
}

/* renumber all exons in a gene */
void FeatureTreePolish::renumberGeneExons(FeatureNode* geneTree) const {
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        renumberTranscriptExons(geneTree->fChildren[i]);
    }
}

/* Is a feature node remapped */
bool FeatureTreePolish::isRemapped(FeatureNode* featureNode) const {
    return (featureNode->fFeature->findAttr(REMAP_STATUS_ATTR) != NULL)
        and (featureNode->fFeature->findAttr(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR) == NULL);
}

/* add a mapping version to an id */
void FeatureTreePolish::setMappingVersionInId(FeatureNode* featureNode,
                                              const AttrVal* attr,
                                              int version) const {
    featureNode->fFeature->getAttrs().update(AttrVal(attr->getName(),
                                                     attr->getVal() + "_" + toString(version)));
}

/* set mapping version in attributes */
void FeatureTreePolish::setMappingVersion(FeatureNode* featureNode,
                                          const string& idAttrName,
                                          const string& havanaIdAttrName,
                                          int version) const {
    const AttrVal* idAttr = featureNode->fFeature->getAttr(idAttrName);
    setMappingVersionInId(featureNode, idAttr, version);
    const AttrVal* havanaIdAttr = featureNode->fFeature->findAttr(havanaIdAttrName);
    if (havanaIdAttr != NULL) {
        setMappingVersionInId(featureNode, havanaIdAttr, version);
    }
}

/* recursively set mapping version in attributes */
void FeatureTreePolish::recursiveSetMappingVersion(FeatureNode* featureNode,
                                                   const string& idAttrName,
                                                   const string& havanaIdAttrName,
                                                   int version) const {
    setMappingVersion(featureNode, idAttrName, havanaIdAttrName, version);
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        recursiveSetMappingVersion(featureNode->fChildren[i], idAttrName, havanaIdAttrName, version);
    }
}

/* Added mapping version numbers */
void FeatureTreePolish::setTranscriptMappingVersion(FeatureNode* transcriptTree) const {
    int version = 1; // FIXME: tmp
    recursiveSetMappingVersion(transcriptTree, GxfFeature::TRANSCRIPT_ID_ATTR, GxfFeature::TRANSCRIPT_HAVANA_ATTR, version);
    // FIXME: exon_id
}

/* Added mapping version numbers */
void FeatureTreePolish::setGeneMappingVersion(FeatureNode* geneTree) const {
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        if (isRemapped(geneTree->fChildren[i])) {
            setTranscriptMappingVersion(geneTree->fChildren[i]);
        }
    }
    int version = 1; // FIXME: tmp
    recursiveSetMappingVersion(geneTree, GxfFeature::GENE_ID_ATTR, GxfFeature::GENE_HAVANA_ATTR, version);
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* geneTree) const {
    if (isRemapped(geneTree)) {
        setGeneMappingVersion(geneTree);
    }
    renumberGeneExons(geneTree);
}


