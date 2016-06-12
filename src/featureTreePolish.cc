/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>
#include <algorithm>
#include "annotationSet.hh"

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
bool FeatureTreePolish::isRemapped(const FeatureNode* featureNode) const {
    return (featureNode->fFeature->findAttr(REMAP_STATUS_ATTR) != NULL)
        and (featureNode->fFeature->findAttr(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR) == NULL);
}

/* find a gene or transcript previously mapped feature by id, or NULL if not
 * mapped or no previous */
const FeatureNode* FeatureTreePolish::getPrevMappedFeature(const FeatureNode* newFeature) const {
    if (fPreviousMappedAnotations == NULL) {
        return NULL;
    } else {
        return fPreviousMappedAnotations->getFeatureNodeById(newFeature->getTypeId(),
                                                             newFeature->fFeature->fSeqid);
    }
}

/* compute mapping version given a possible-null previous feature and
 * if it's considered the same */
int FeatureTreePolish::getFeatureMappingVersion(const FeatureNode* prevFeature,
                                                bool featureSame) const {
    if (prevFeature == NULL) {
        return 1;
    } else {
        int prevMappingVersion = getMappingVersion(prevFeature->getTypeId());
        // if there is no previous mapping version, which either means the
        // set was created before we had mapping versions or the previous
        // version was target substituted.
        // FIXME: really need to keep complete version history as an
        // entry could switch between mapped and substituted
        if (prevMappingVersion == 0) {
            return featureSame ? 1 : 2;
        } else {
            return featureSame ? prevMappingVersion : prevMappingVersion+1;
        }
    }
}

/* compare two node attribute values */
bool FeatureTreePolish::compareNodeAttrVals(const FeatureNode* prevNode,
                                            const FeatureNode* newNode,
                                            const string& attrName) const {
    const AttrVal* prevAttr = prevNode->fFeature->getAttrs().find(attrName);
    const AttrVal* newAttr = newNode->fFeature->getAttrs().find(attrName);
    if ((prevAttr == NULL) and (newAttr == NULL)) {
        return true;  // both NULL
    } else if ((prevAttr == NULL) or (newAttr == NULL)) {
        return false;  // one NULL
    } else if (prevAttr->getVals().size() != newAttr->getVals().size()) {
        return false;  // different number of values
    } else if (prevAttr->getVals().size() == 1) {
        return (prevAttr->getVal() == newAttr->getVal());  // fast path
    } else {
        // sort ignore order
        StringVector prevVals(prevAttr->getVals());
        sort(prevVals.begin(), prevVals.end());
        StringVector newVals(newAttr->getVals());
        sort(newVals.begin(), newVals.end());
        return (prevVals == newVals);
    }
}

/* compare a list of two node attribute/values */
bool FeatureTreePolish::compareNodeAttrs(const FeatureNode* prevNode,
                                         const FeatureNode* newNode,
                                         const StringVector& attrNames) const {
    for (int i = 0; i < attrNames.size(); i++) {
        if (not compareNodeAttrVals(prevNode, newNode, attrNames[i])) {
            return false;
        }
    }
    return true;
}


/* compare a mapped node with previous mapped node.  This is not recursive */
bool FeatureTreePolish::compareMappedNodes(const FeatureNode* prevNode,
                                           const FeatureNode* newNode,
                                           const StringVector& attrNames) const {
    const GxfFeature* prevFeat = prevNode->fFeature;
    const GxfFeature* newFeat = newNode->fFeature;
    return (prevFeat->fSource == newFeat->fSource)
        and (prevFeat->fStart == newFeat->fStart)
        and (prevFeat->fEnd == newFeat->fEnd)
        and (prevFeat->fStrand == newFeat->fStrand)
        and (prevFeat->fPhase == newFeat->fPhase)
        and compareNodeAttrs(prevNode, newNode, attrNames);
}

/* compare a mapped gene node with previous mapped gene.  This is not
 * recursive */
bool FeatureTreePolish::compareGeneNodes(const FeatureNode* prevNode,
                                         const FeatureNode* newNode) const {
    static const StringVector attrNames = {
        GxfFeature::GENE_NAME_ATTR,
        GxfFeature::GENE_TYPE_ATTR,
        GxfFeature::GENE_STATUS_ATTR,
        GxfFeature::GENE_HAVANA_ATTR,
        GxfFeature::TAG_ATTR
    };
    return compareMappedNodes(prevNode, newNode, attrNames);
}

/* compare a mapped transcript node with previous mapped transcript.  This is not
 * recursive */
bool FeatureTreePolish::compareTranscriptNodes(const FeatureNode* prevNode,
                                               const FeatureNode* newNode) const {
    static const StringVector attrNames = {
        GxfFeature::GENE_NAME_ATTR,
        GxfFeature::TRANSCRIPT_NAME_ATTR,
        GxfFeature::TRANSCRIPT_TYPE_ATTR,
        GxfFeature::TRANSCRIPT_STATUS_ATTR,
        GxfFeature::TRANSCRIPT_HAVANA_ATTR,
        GxfFeature::TAG_ATTR
    };
    return compareMappedNodes(prevNode, newNode, attrNames);
}

/* compare a mapped node of other with previous mapped node.  This is not
 * recursive */
bool FeatureTreePolish::compareOtherNodes(const FeatureNode* prevNode,
                                          const FeatureNode* newNode) const {
    static const StringVector attrNames = {
        GxfFeature::EXON_ID_ATTR,
        GxfFeature::EXON_NUMBER_ATTR,
        GxfFeature::TAG_ATTR
    };
    return compareMappedNodes(prevNode, newNode, attrNames);
}

/* recursively compare descendant nodes of a transcript with previous mapped
 * transcript.  Parent should have already been checked */
bool FeatureTreePolish::compareMappedTranscriptsDescendants(const FeatureNode* prevParent,
                                                            const FeatureNode* newParent) const {
    if (prevParent->fChildren.size() != newParent->fChildren.size()) {
        return false;
    } else if (prevParent->fChildren.size() == 0) {
        return true;
    } else {
        // compare children sorted 
        FeatureNodeVector prevChildren(prevParent->fChildren);
        prevChildren.sort();
        FeatureNodeVector newChildren(newParent->fChildren);
        newChildren.sort();

        for (int i = 0; i < prevChildren.size(); i++) {
            if (not compareOtherNodes(prevChildren[i], newChildren[i])) {
                return false;
            }
            if (not compareMappedTranscriptsDescendants(prevChildren[i], newChildren[i])) {
                return false;
            }
        }
        return true;
    }
}

/* recursively compare nodes of a transcript with previous mapped transcript. */
bool FeatureTreePolish::compareMappedTranscripts(const FeatureNode* prevTranscript,
                                                 const FeatureNode* newTranscript) const {
    assert(prevTranscript->fFeature->fType == GxfFeature::TRANSCRIPT);

    if (not compareTranscriptNodes(prevTranscript, newTranscript)) {
        return false;
    } else {
        return compareMappedTranscriptsDescendants(prevTranscript, newTranscript);
    }
}

/* add a mapping version to an id */
void FeatureTreePolish::setMappingVersionInId(FeatureNode* featureNode,
                                              const AttrVal* attr,
                                              int mappingVersion) const {
    featureNode->fFeature->getAttrs().update(AttrVal(attr->getName(),
                                                     attr->getVal() + "_" + toString(mappingVersion)));
}

/* set mapping version in attributes */
void FeatureTreePolish::setMappingVersion(FeatureNode* featureNode,
                                          const string& idAttrName,
                                          const string& havanaIdAttrName,
                                          int mappingVersion) const {
    const AttrVal* idAttr = featureNode->fFeature->getAttr(idAttrName);
    setMappingVersionInId(featureNode, idAttr, mappingVersion);
    const AttrVal* havanaIdAttr = featureNode->fFeature->findAttr(havanaIdAttrName);
    if (havanaIdAttr != NULL) {
        setMappingVersionInId(featureNode, havanaIdAttr, mappingVersion);
    }
}

/* recursively set mapping version in attributes */
void FeatureTreePolish::recursiveSetMappingVersion(FeatureNode* featureNode,
                                                   const string& idAttrName,
                                                   const string& havanaIdAttrName,
                                                   int mappingVersion) const {
    setMappingVersion(featureNode, idAttrName, havanaIdAttrName, mappingVersion);
    for (int i = 0; i < featureNode->fChildren.size(); i++) {
        recursiveSetMappingVersion(featureNode->fChildren[i], idAttrName, havanaIdAttrName, mappingVersion);
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

/* Added mapping version numbers.  Return true if transcript is the same
 * as the previous or new, or false if it has changed */
bool FeatureTreePolish::setTranscriptMappingVersion(FeatureNode* transcriptTree) const {
    // find previous transcript, if it exists and derive version from it.
    const FeatureNode* prevTranscript = getPrevMappedFeature(transcriptTree);
    bool transcriptSame = (prevTranscript == NULL)
        or compareMappedTranscriptsDescendants(transcriptTree, prevTranscript);
    int mappingVersion = getFeatureMappingVersion(prevTranscript, transcriptSame);
  
    recursiveSetMappingVersion(transcriptTree, GxfFeature::TRANSCRIPT_ID_ATTR, GxfFeature::TRANSCRIPT_HAVANA_ATTR, mappingVersion);
    return transcriptSame;
}

/* Added mapping version numbers to a exon feature */
void FeatureTreePolish::setExonMappingVersion(FeatureNode* exonFeature,
                                              int mappingVersion) const {
    const AttrVal* idAttr = exonFeature->fFeature->getAttr(GxfFeature::EXON_ID_ATTR);
    setMappingVersionInId(exonFeature, idAttr, mappingVersion);
}

/* Added mapping version numbers to a set of exons feeatures with the same id */
void FeatureTreePolish::setExonMappingVersion(const string& exonId,
                                              vector<FeatureNode*> exonFeatures) const {
    int mappingVersion = 1; // FIXME: tmp
    for (int i = 0; i < exonFeatures.size(); i++) {
        setExonMappingVersion(exonFeatures[i], mappingVersion);
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
            setTranscriptMappingVersion(geneTree->fChildren[i]);
            recordTranscriptMappedExons(geneTree->fChildren[i], exonIdExonMap);
        }
    }
    int mappingVersion = 1; // FIXME: tmp
    recursiveSetMappingVersion(geneTree, GxfFeature::GENE_ID_ATTR, GxfFeature::GENE_HAVANA_ATTR, mappingVersion);
    setExonsMappingVersions(exonIdExonMap);
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* geneTree) const {
    if (isRemapped(geneTree)) {
        setGeneMappingVersion(geneTree);
    }
    renumberGeneExons(geneTree);
}


