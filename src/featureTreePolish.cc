/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>
#include <algorithm>
#include "annotationSet.hh"
#include <iostream>

const bool DEBUG = true;

/* renumber ab exon */
void FeatureTreePolish::renumberExon(Feature* exon,
                                     int exonNum,
                                     ExonNumExonMap& exonNumExonMap) const {
    int oldExonNum = stringToInt(exon->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR)->getVal());
    exonNumExonMap[oldExonNum].push_back(exon);
    exon->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
void FeatureTreePolish::renumberExons(Feature* transcript,
                                      ExonNumExonMap& exonNumExonMap) const {
    int exonNum = 1;
    // exons are always in genomic order.
    for (int i = 0; i < transcript->getChildren().size(); i++) {
        if (transcript->getChild(i)->getType() == GxfFeature::EXON) {
            renumberExon(transcript->getChild(i), exonNum++, exonNumExonMap);
        }
    }
}

/* find the new exon containing a feature given the old exon number */
Feature* FeatureTreePolish::findNewExon(Feature* feature,
                                        int oldExonNum,
                                        ExonNumExonMap& exonNumExonMap) const {
    for (int i = 0; i < exonNumExonMap[oldExonNum].size(); i++) {
        if (feature->overlaps(exonNumExonMap[oldExonNum][i])) {
            return exonNumExonMap[oldExonNum][i];
        }
    }
    throw logic_error("renumberOtherFeature: lost exon");
}

/* change a non-exon feature to their exon numbers */
void FeatureTreePolish::renumberOtherFeature(Feature* feature,
                                             ExonNumExonMap& exonNumExonMap) const {
    const AttrVal* exonNumAttr = feature->getAttrs().find(GxfFeature::EXON_NUMBER_ATTR);;
    if (exonNumAttr != NULL) {
        Feature *newExon = findNewExon(feature, stringToInt(exonNumAttr->getVal()), exonNumExonMap);
        const AttrVal* newExonNumAttr = newExon->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR);
        feature->getAttrs().update(*newExonNumAttr);
    }
}

/* recursively change non-exons features to match the news exon numbers */
void FeatureTreePolish::renumberOtherFeatures(Feature* feature,
                                              ExonNumExonMap& exonNumExonMap) const {
    for (int i = 0; i < feature->getChildren().size(); i++) {
        if (feature->getChild(i)->getType() != GxfFeature::EXON) {
            renumberOtherFeature(feature->getChild(i), exonNumExonMap);
        }
        renumberOtherFeatures(feature->getChild(i), exonNumExonMap);
    }
}

/* renumber all features in a transcript */
void FeatureTreePolish::renumberTranscriptExons(Feature* transcript) const {
    assert(transcript->getType() == GxfFeature::TRANSCRIPT);
    ExonNumExonMap exonNumExonMap;
    renumberExons(transcript, exonNumExonMap);
    renumberOtherFeatures(transcript, exonNumExonMap);
}

/* renumber all exons in a gene */
void FeatureTreePolish::renumberGeneExons(Feature* gene) const {
    for (int i = 0; i < gene->getChildren().size(); i++) {
        renumberTranscriptExons(gene->getChild(i));
    }
}

/* Is a feature remapped */
bool FeatureTreePolish::isRemapped(const Feature* feature) const {
    return (feature->findAttr(REMAP_STATUS_ATTR) != NULL)
        and (feature->findAttr(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR) == NULL);
}

/* find a gene or transcript previously mapped feature by id, or NULL if not
 * mapped or no previous */
const Feature* FeatureTreePolish::getPrevMappedFeature(const Feature* newFeature) const {
    if (fPreviousMappedAnotations == NULL) {
        return NULL;
    } else {
        return fPreviousMappedAnotations->getFeatureById(newFeature->getTypeId(),
                                                         newFeature->getSeqid());
    }
}

/* compute mapping version given a possible-null previous feature and
 * if it's considered the same */
int FeatureTreePolish::getFeatureMappingVersion(const Feature* prevFeature,
                                                bool featureSame) const {
    if ((prevFeature == NULL) || (not isRemapped(prevFeature))) {
        return 1;  // first mapped version
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

/* compare one value */
bool FeatureTreePolish::compareAttrVal(const AttrVal* prevAttr,
                                       const AttrVal* newAttr,
                                       int iValue,
                                       bool isIdAttr) const {
    if (isIdAttr) {
        return (getPreMappedId(prevAttr->getVal(iValue)) == getPreMappedId(newAttr->getVal(iValue)));
    } else {
        return (prevAttr->getVal(iValue) == newAttr->getVal(iValue));
    }
}


/* compare two feature attribute values. If it's an id attribute, it's compared 
* without mapping versions */
bool FeatureTreePolish::compareAttrVals(const Feature* prevFeature,
                                        const Feature* newFeature,
                                        const string& attrName,
                                        bool isIdAttr) const {
    const AttrVal* prevAttr = prevFeature->getAttrs().find(attrName);
    const AttrVal* newAttr = newFeature->getAttrs().find(attrName);
    if ((prevAttr == NULL) and (newAttr == NULL)) {
        return true;  // both NULL
    } else if ((prevAttr == NULL) or (newAttr == NULL)) {
        return false;  // one NULL
    } else if (prevAttr->getVals().size() != newAttr->getVals().size()) {
        return false;  // different number of values
    } else {
        for (int iVal = 0; iVal < prevAttr->getVals().size(); iVal++) {
            if (not compareAttrVal(prevAttr, newAttr, iVal, isIdAttr)) {
                return false;
            }
        }
        return true;
    }
}

/* compare a list of two feature attribute/values */
bool FeatureTreePolish::compareAttrs(const Feature* prevFeature,
                                     const Feature* newFeature,
                                     const StringVector& attrNames,
                                     bool isIdAttrs) const {
    for (int i = 0; i < attrNames.size(); i++) {
        if (not compareAttrVals(prevFeature, newFeature, attrNames[i], isIdAttrs)) {
            return false;
        }
    }
    return true;
}

/* compare a list of two feature attribute/values */
bool FeatureTreePolish::compareAttrs(const Feature* prevFeature,
                                     const Feature* newFeature,
                                     const StringVector& attrNames,
                                     const StringVector& idAttrNames) const {
    return compareAttrs(prevFeature, newFeature, attrNames, false)
        and compareAttrs(prevFeature, newFeature, idAttrNames, true);
}


/* compare a mapped node with previous mapped node.  This is not recursive */
bool FeatureTreePolish::compareMappedFeatures(const Feature* prevFeature,
                                              const Feature* newFeature,
                                              const StringVector& attrNames,
                                              const StringVector& idAttrNames) const {
    bool same = (prevFeature->getSource() == newFeature->getSource())
        and (prevFeature->getStart() == newFeature->getStart())
        and (prevFeature->getEnd() == newFeature->getEnd())
        and (prevFeature->getStrand() == newFeature->getStrand())
        and (prevFeature->getPhase() == newFeature->getPhase())
        and compareAttrs(prevFeature, newFeature, attrNames, idAttrNames);
    if (DEBUG and not same) {
        cerr << "diff prev\t" << prevFeature->toString() << endl
             << "      new\t" << newFeature->toString() << endl;
    }
    return same;
}

/* compare a mapped gene node with previous mapped gene.  This is not
 * recursive */
bool FeatureTreePolish::compareGeneFeatures(const Feature* prevFeature,
                                            const Feature* newFeature) const {
    static const StringVector attrNames = {
        GxfFeature::GENE_NAME_ATTR,
        GxfFeature::GENE_TYPE_ATTR,
        GxfFeature::GENE_STATUS_ATTR,
        GxfFeature::TAG_ATTR
    };
    static const StringVector idAttrNames = {
        GxfFeature::GENE_ID_ATTR,
        GxfFeature::GENE_HAVANA_ATTR,
    };
    return compareMappedFeatures(prevFeature, newFeature, attrNames, idAttrNames);
}

/* compare a mapped transcript node with previous mapped transcript.  This is not
 * recursive */
bool FeatureTreePolish::compareTranscriptFeatures(const Feature* prevFeature,
                                                  const Feature* newFeature) const {
    static const StringVector attrNames = {
        GxfFeature::GENE_NAME_ATTR,
        GxfFeature::TRANSCRIPT_NAME_ATTR,
        GxfFeature::TRANSCRIPT_TYPE_ATTR,
        GxfFeature::TRANSCRIPT_STATUS_ATTR,
        GxfFeature::TRANSCRIPT_HAVANA_ATTR,
        GxfFeature::TAG_ATTR
    };
    static const StringVector idAttrNames = {
        GxfFeature::GENE_ID_ATTR,
        GxfFeature::TRANSCRIPT_HAVANA_ATTR,
    };
    return compareMappedFeatures(prevFeature, newFeature, attrNames, idAttrNames);
}

/* compare a mapped node of other with previous mapped node.  This is not
 * recursive */
bool FeatureTreePolish::compareOtherFeatures(const Feature* prevFeature,
                                             const Feature* newFeature) const {
    static const StringVector attrNames = {
        GxfFeature::EXON_NUMBER_ATTR,
        GxfFeature::TAG_ATTR
    };
    static const StringVector idAttrNames = {
        GxfFeature::EXON_ID_ATTR,
    };
    return compareMappedFeatures(prevFeature, newFeature, attrNames, idAttrNames);
}

/* recursively compare descendant nodes of a transcript with previous mapped
 * transcript.  Parent should have already been checked */
bool FeatureTreePolish::compareMappedTranscriptsDescendants(const Feature* prevParent,
                                                            const Feature* newParent) const {
    if (prevParent->getChildren().size() != newParent->getChildren().size()) {
        return false;
    } else if (prevParent->getChildren().size() == 0) {
        return true;
    } else {
        // compare children sorted 
        FeatureVector prevChildren(prevParent->getChildren());
        prevChildren.sort();
        FeatureVector newChildren(newParent->getChildren());
        newChildren.sort();

        for (int i = 0; i < prevChildren.size(); i++) {
            if (not compareOtherFeatures(prevChildren[i], newChildren[i])) {
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
bool FeatureTreePolish::compareMappedTranscripts(const Feature* prevTranscript,
                                                 const Feature* newTranscript) const {
    assert(prevTranscript->getType() == GxfFeature::TRANSCRIPT);

    if (not compareTranscriptFeatures(prevTranscript, newTranscript)) {
        return false;
    } else {
        return compareMappedTranscriptsDescendants(prevTranscript, newTranscript);
    }
}

/* add a mapping version to an id */
void FeatureTreePolish::setMappingVersionInId(Feature* feature,
                                              const AttrVal* attr,
                                              int mappingVersion) const {
    feature->getAttrs().update(AttrVal(attr->getName(),
                                       attr->getVal() + "_" + toString(mappingVersion)));
}

/* set mapping version in attributes */
void FeatureTreePolish::setMappingVersion(Feature* feature,
                                          const string& idAttrName,
                                          const string& havanaIdAttrName,
                                          int mappingVersion) const {
    const AttrVal* idAttr = feature->getAttr(idAttrName);
    setMappingVersionInId(feature, idAttr, mappingVersion);
    const AttrVal* havanaIdAttr = feature->findAttr(havanaIdAttrName);
    if (havanaIdAttr != NULL) {
        setMappingVersionInId(feature, havanaIdAttr, mappingVersion);
    }
}

/* recursively set mapping version in attributes */
void FeatureTreePolish::recursiveSetMappingVersion(Feature* feature,
                                                   const string& idAttrName,
                                                   const string& havanaIdAttrName,
                                                   int mappingVersion) const {
    setMappingVersion(feature, idAttrName, havanaIdAttrName, mappingVersion);
    for (int i = 0; i < feature->getChildren().size(); i++) {
        recursiveSetMappingVersion(feature->getChild(i), idAttrName, havanaIdAttrName, mappingVersion);
    }
}
/*
 * record remapped exons by id for latter adding mapping version.
 */
void FeatureTreePolish::recordTranscriptMappedExons(Feature* transcript,
                                                    ExonIdExonMap& exonIdExonMap) const {
    for (int i = 0; i < transcript->getChildren().size(); i++) {
        Feature* child = transcript->getChild(i);
        if (child->isExon()) {
            exonIdExonMap[child->getTypeId()].push_back(child);
        }
    }
}

/* Added mapping version numbers.  Return true if transcript is the same
 * as the previous or new, or false if it has changed */
bool FeatureTreePolish::setTranscriptMappingVersion(Feature* transcript) const {
    // find previous transcript, if it exists and derive version from it.
    const Feature* prevTranscript = getPrevMappedFeature(transcript);
    bool transcriptSame = (prevTranscript == NULL)
        or compareMappedTranscriptsDescendants(prevTranscript, transcript);
    int mappingVersion = getFeatureMappingVersion(prevTranscript, transcriptSame);
  
    recursiveSetMappingVersion(transcript, GxfFeature::TRANSCRIPT_ID_ATTR, GxfFeature::TRANSCRIPT_HAVANA_ATTR, mappingVersion);
    return transcriptSame;
}

/* Added mapping version numbers to a exon feature */
void FeatureTreePolish::setExonMappingVersion(Feature* exonFeature,
                                              int mappingVersion) const {
    const AttrVal* idAttr = exonFeature->getAttr(GxfFeature::EXON_ID_ATTR);
    setMappingVersionInId(exonFeature, idAttr, mappingVersion);
}

/* Added mapping version numbers to a set of exons feeatures with the same id */
void FeatureTreePolish::setExonMappingVersion(const string& exonId,
                                              vector<Feature*> exonFeatures) const {
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
void FeatureTreePolish::setGeneMappingVersion(Feature* gene) const {
    // exon id scope is gene
    ExonIdExonMap exonIdExonMap;
    for (int i = 0; i < gene->getChildren().size(); i++) {
        if (isRemapped(gene->getChild(i))) {
            setTranscriptMappingVersion(gene->getChild(i));
            recordTranscriptMappedExons(gene->getChild(i), exonIdExonMap);
        }
    }
    int mappingVersion = 1; // FIXME: tmp
    recursiveSetMappingVersion(gene, GxfFeature::GENE_ID_ATTR, GxfFeature::GENE_HAVANA_ATTR, mappingVersion);
    setExonsMappingVersions(exonIdExonMap);
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(Feature* gene) const {
    // n.b. must renumber exons first, otherwise different exon numbers
    // with previous version will cause false mapping version increments.
    renumberGeneExons(gene);
    if (isRemapped(gene)) {
        setGeneMappingVersion(gene);
    }
}


