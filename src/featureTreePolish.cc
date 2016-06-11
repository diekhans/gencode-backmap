/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>

/* renumber ab exon */
void FeatureTreePolish::renumberExon(FeatureNode* exonNode,
                                     int exonNum,
                                     ExonNumExonMap& exonNumExonMap) const {
    int oldExonNum = stringToInt(exonNode->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR)->getVal());
    exonNumExonMap[oldExonNum].push_back(exonNode);
    exonNode->fFeature->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
void FeatureTreePolish::renumberExons(FeatureNode* transcriptTree,
                                      ExonNumExonMap& exonNumExonMap) const {
    int exonNum = 1;
    // exons are always in genomic order.
    for (int i = 0; i < transcriptTree->fChildren.size(); i++) {
        if (transcriptTree->fChildren[i]->fFeature->fType == GxfFeature::EXON) {
            renumberExon(transcriptTree->fChildren[i], exonNum++, exonNumExonMap);
        }
    }
}

/* find the new exon containing a feature given the old exon number */
FeatureNode* FeatureTreePolish::findNewExon(FeatureNode* featureNode,
                                            int oldExonNum,
                                            ExonNumExonMap& exonNumExonMap) const {
    for (int i = 0; i < exonNumExonMap[oldExonNum].size(); i++) {
        if (featureNode->fFeature->overlaps(exonNumExonMap[oldExonNum][i]->fFeature)) {
            return exonNumExonMap[oldExonNum][i];
        }
    }
    throw logic_error("renumberOtherFeature: lost exon");
}

/* change a non-exon feature to their exon numbers */
void FeatureTreePolish::renumberOtherFeature(FeatureNode* featureNode,
                                             ExonNumExonMap& exonNumExonMap) const {
    const AttrVal* exonNumAttr = featureNode->fFeature->getAttrs().find(GxfFeature::EXON_NUMBER_ATTR);;
    if (exonNumAttr != NULL) {
        FeatureNode *newExon = findNewExon(featureNode, stringToInt(exonNumAttr->getVal()), exonNumExonMap);
        const AttrVal* newExonNumAttr = newExon->fFeature->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR);
        featureNode->fFeature->getAttrs().update(*newExonNumAttr);
    }
}

/* recursively change non-exons features to match the news exon numbers */
void FeatureTreePolish::renumberOtherFeatures(FeatureNode* featureNode,
                                              ExonNumExonMap& exonNumExonMap) const {
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        if (featureNode->fChildren[i]->fFeature->fType != GxfFeature::EXON) {
            renumberOtherFeature(featureNode->fChildren[i], exonNumExonMap);
        }
        renumberOtherFeatures(featureNode->fChildren[i], exonNumExonMap);
    }
}

/* renumber all features in a transcript */
void FeatureTreePolish::renumberTranscriptExons(FeatureNode* transcriptTree) const {
    assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);
    ExonNumExonMap exonNumExonMap;
    renumberExons(transcriptTree, exonNumExonMap);
    renumberOtherFeatures(transcriptTree, exonNumExonMap);
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

/*
 * record remapped exons by id for latter adding mapping version.
 */
void FeatureTreePolish::recordTranscriptMappedExons(FeatureNode* transcriptTree,
                                                    ExonIdExonMap& exonIdExonMap) const {
    for (int i = 0; i < transcriptTree->fChildren.size(); i++) {
        FeatureNode* child = transcriptTree->fChildren[i];
        if (child->isExon()) {
            exonIdExonMap[child->getTypeId()].push_back(child);
        }
    }
}

/* Added mapping version numbers */
void FeatureTreePolish::setTranscriptMappingVersion(FeatureNode* transcriptTree,
                                                    ExonIdExonMap& exonIdExonMap) const {
    int version = 1; // FIXME: tmp
    recursiveSetMappingVersion(transcriptTree, GxfFeature::TRANSCRIPT_ID_ATTR, GxfFeature::TRANSCRIPT_HAVANA_ATTR, version);
    recordTranscriptMappedExons(transcriptTree, exonIdExonMap);
    // FIXME: exon_id
}

/* Added mapping version numbers to a exon feature */
void FeatureTreePolish::setExonMappingVersion(FeatureNode* exonFeature,
                                              int version) const {
    const AttrVal* idAttr = exonFeature->fFeature->getAttr(GxfFeature::EXON_ID_ATTR);
    setMappingVersionInId(exonFeature, idAttr, version);
}

/* Added mapping version numbers to a set of exons feeatures with the same id */
void FeatureTreePolish::setExonMappingVersion(const string& exonId,
                                              vector<FeatureNode*> exonFeatures) const {
    int version = 1; // FIXME: tmp
    for (int i = 0; i < exonFeatures.size(); i++) {
        setExonMappingVersion(exonFeatures[i], version);
    }
}

/* Added mapping version numbers to all exons */
void FeatureTreePolish::setExonsMappingVersions(ExonIdExonMap& exonIdExonMap) const {
    for (ExonIdExonMapIter iter = exonIdExonMap.begin(); iter != exonIdExonMap.end(); iter++) {
        setExonMappingVersion(iter->first, iter->second);
    }
}

/* Added mapping version numbers */
void FeatureTreePolish::setGeneMappingVersion(FeatureNode* geneTree) const {
    // exon id is scope is-gene
    ExonIdExonMap exonIdExonMap;
    for (int i = 0; i < geneTree->fChildren.size(); i++) {
        if (isRemapped(geneTree->fChildren[i])) {
            setTranscriptMappingVersion(geneTree->fChildren[i], exonIdExonMap);
        }
    }
    int version = 1; // FIXME: tmp
    recursiveSetMappingVersion(geneTree, GxfFeature::GENE_ID_ATTR, GxfFeature::GENE_HAVANA_ATTR, version);
    setExonsMappingVersions(exonIdExonMap);
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* geneTree) const {
    if (isRemapped(geneTree)) {
        setGeneMappingVersion(geneTree);
    }
    renumberGeneExons(geneTree);
}


