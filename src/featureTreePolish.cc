/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>
#include <algorithm>
#include "annotationSet.hh"
#include <iostream>
#include "frame.hh"

static const bool DEBUG = false;

/* renumber an exon */
static void renumberExon(FeatureNode* exon,
                         int exonNum,
                         FeatureTreePolish::ExonNumExonMap& exonNumExonMap) {
    int oldExonNum = stringToInt(exon->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR)->getVal());
    exonNumExonMap[oldExonNum].push_back(exon);
    exon->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* renumber all exons in a transcript */
static void renumberExons(FeatureNode* transcript,
                          FeatureTreePolish::ExonNumExonMap& exonNumExonMap) {
    int exonNum = 1;
    // exons are always in genomic order.
    for (int i = 0; i < transcript->getNumChildren(); i++) {
        if (transcript->getChild(i)->getType() == GxfFeature::EXON) {
            renumberExon(transcript->getChild(i), exonNum++, exonNumExonMap);
        }
    }
}

/* find the new exon containing a feature given the old exon number */
static FeatureNode* findNewExon(FeatureNode* feature,
                                int oldExonNum,
                                FeatureTreePolish::ExonNumExonMap& exonNumExonMap) {
    for (int i = 0; i < exonNumExonMap[oldExonNum].size(); i++) {
        if (feature->overlaps(exonNumExonMap[oldExonNum][i])) {
            return exonNumExonMap[oldExonNum][i];
        }
    }
    throw logic_error("renumberOtherFeature: lost exon");
}

/* change a non-exon feature to their exon numbers */
static void renumberOtherFeature(FeatureNode* feature,
                                 FeatureTreePolish::ExonNumExonMap& exonNumExonMap) {
    const AttrVal* exonNumAttr = feature->getAttrs().find(GxfFeature::EXON_NUMBER_ATTR);;
    if (exonNumAttr != NULL) {
        FeatureNode* newExon = findNewExon(feature, stringToInt(exonNumAttr->getVal()), exonNumExonMap);
        const AttrVal* newExonNumAttr = newExon->getAttrs().get(GxfFeature::EXON_NUMBER_ATTR);
        feature->getAttrs().update(*newExonNumAttr);
    }
}

/* recursively change non-exons features to match the news exon numbers */
static void renumberOtherFeatures(FeatureNode* feature,
                                  FeatureTreePolish::ExonNumExonMap& exonNumExonMap) {
    for (int i = 0; i < feature->getNumChildren(); i++) {
        if (feature->getChild(i)->getType() != GxfFeature::EXON) {
            renumberOtherFeature(feature->getChild(i), exonNumExonMap);
        }
        renumberOtherFeatures(feature->getChild(i), exonNumExonMap);
    }
}

/* renumber all features in a transcript */
static void renumberTranscriptExons(FeatureNode* transcript) {
    assert(transcript->getType() == GxfFeature::TRANSCRIPT);
    FeatureTreePolish::ExonNumExonMap exonNumExonMap;
    renumberExons(transcript, exonNumExonMap);
    renumberOtherFeatures(transcript, exonNumExonMap);
}

/* renumber all exons in a gene */
static void renumberGeneExons(FeatureNode* gene) {
    for (int i = 0; i < gene->getNumChildren(); i++) {
        renumberTranscriptExons(gene->getChild(i));
    }
}

/* Is a feature remapped */
static bool isRemapped(const FeatureNode* feature) {
    return (feature->findAttr(REMAP_STATUS_ATTR) != NULL)
        and (feature->findAttr(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR) == NULL);
}

/* compute mapping version given a possible-null previous feature and
 * if it's considered the same */
static int getFeatureMappingVersion(const FeatureNode* prevFeature,
                                    bool featureSame) {
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
static bool compareAttrVal(const AttrVal* prevAttr,
                           const AttrVal* newAttr,
                           int iValue,
                           bool isIdAttr) {
    if (isIdAttr) {
        return (getPreMappedId(prevAttr->getVal(iValue)) == getPreMappedId(newAttr->getVal(iValue)));
    } else {
        return (prevAttr->getVal(iValue) == newAttr->getVal(iValue));
    }
}


/* compare two feature attribute values. If it's an id attribute, it's compared 
* without mapping versions */
static bool compareAttrVals(const FeatureNode* prevFeature,
                            const FeatureNode* newFeature,
                            const string& attrName,
                            bool isIdAttr) {
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
static bool compareAttrs(const FeatureNode* prevFeature,
                         const FeatureNode* newFeature,
                         const StringVector& attrNames,
                         bool isIdAttrs) {
    for (int i = 0; i < attrNames.size(); i++) {
        if (not compareAttrVals(prevFeature, newFeature, attrNames[i], isIdAttrs)) {
            return false;
        }
    }
    return true;
}

/* compare a list of two feature attribute/values */
static bool compareAttrs(const FeatureNode* prevFeature,
                         const FeatureNode* newFeature,
                         const StringVector& attrNames,
                         const StringVector& idAttrNames) {
    return compareAttrs(prevFeature, newFeature, attrNames, false)
        and compareAttrs(prevFeature, newFeature, idAttrNames, true);
}


/* compare a mapped feature with previous mapped feature.  This is not recursive */
static bool compareMappedFeatures(const FeatureNode* prevFeature,
                                  const FeatureNode* newFeature,
                                  const StringVector& attrNames,
                                  const StringVector& idAttrNames) {
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

/* compare a mapped gene Feature with previous mapped gene.  This is not
 * recursive */
static bool compareGeneFeatures(const FeatureNode* prevFeature,
                                const FeatureNode* newFeature) {
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

/* compare a mapped feature of other types (CDS, exon, etc). with previous
 * mapped features.  This is not recursive */
static bool compareOtherFeatures(const FeatureNode* prevFeature,
                                             const FeatureNode* newFeature) {
    static const StringVector attrNames = {
        GxfFeature::EXON_NUMBER_ATTR,
        GxfFeature::TAG_ATTR
    };
    static const StringVector idAttrNames = {
        GxfFeature::EXON_ID_ATTR,
    };
    return compareMappedFeatures(prevFeature, newFeature, attrNames, idAttrNames);
}

/* recursively compare descendant features of a transcript with previous mapped
 * transcript.  Parent should have already been checked */
static bool compareMappedTranscriptsDescendants(const FeatureNode* prevParent,
                                                const FeatureNode* newParent) {
    if (prevParent->getNumChildren() != newParent->getNumChildren()) {
        return false;
    } else if (prevParent->getNumChildren() == 0) {
        return true;
    } else {
        // compare children sorted 
        FeatureNodeVector prevChildren(prevParent->getChildren());
        prevChildren.sort();
        FeatureNodeVector newChildren(newParent->getChildren());
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

/* add a mapping version to an id */
static void setMappingVersionInId(FeatureNode* feature,
                                  const AttrVal* attr,
                                  int mappingVersion) {
    feature->getAttrs().update(AttrVal(attr->getName(),
                                       attr->getVal() + "_" + toString(mappingVersion)));
}

/* set mapping version in attributes */
static void setMappingVersion(FeatureNode* feature,
                              const string& idAttrName,
                              const string& havanaIdAttrName,
                              int mappingVersion) {
    const AttrVal* idAttr = feature->getAttr(idAttrName);
    setMappingVersionInId(feature, idAttr, mappingVersion);
    const AttrVal* havanaIdAttr = feature->findAttr(havanaIdAttrName);
    if (havanaIdAttr != NULL) {
        setMappingVersionInId(feature, havanaIdAttr, mappingVersion);
    }
}

/* recursively set mapping version in attributes */
static void recursiveSetMappingVersion(FeatureNode* feature,
                                       const string& idAttrName,
                                       const string& havanaIdAttrName,
                                       int mappingVersion) {
    setMappingVersion(feature, idAttrName, havanaIdAttrName, mappingVersion);
    for (int i = 0; i < feature->getNumChildren(); i++) {
        recursiveSetMappingVersion(feature->getChild(i), idAttrName, havanaIdAttrName, mappingVersion);
    }
}

/* find a gene or transcript previously mapped feature by id, or NULL if not
 * mapped or no previous */
const FeatureNode* FeatureTreePolish::getPrevMappedFeature(const FeatureNode* newFeature) const {
    if (fPreviousMappedAnotations == NULL) {
        return NULL;
    } else {
        return fPreviousMappedAnotations->getFeatureById(newFeature->getTypeId(),
                                                         newFeature->isParY());
    }
}

/* Added mapping version numbers.  Return true if transcript is the same
 * as the previous or new, or false if it has changed */
bool FeatureTreePolish::setTranscriptMappingVersion(FeatureNode* transcript) const {
    // find previous transcript, if it exists and derive version from it.
    const FeatureNode* prevTranscript = getPrevMappedFeature(transcript);
    bool transcriptSame = (prevTranscript == NULL)
        or compareMappedTranscriptsDescendants(prevTranscript, transcript);
    int mappingVersion = getFeatureMappingVersion(prevTranscript, transcriptSame);
  
    recursiveSetMappingVersion(transcript, GxfFeature::TRANSCRIPT_ID_ATTR, GxfFeature::TRANSCRIPT_HAVANA_ATTR, mappingVersion);
    return transcriptSame;
}

/*
 * recursive collect exons by pre-mapped id.
 */
static void collectExons(const FeatureNode* root,
                         FeatureTreePolish::ExonIdExonMap& exonIdExonMap) {
    for (int i = 0; i < root->getNumChildren(); i++) {
        const FeatureNode* child = root->getChild(i);
        if (child->isExon()) {
            exonIdExonMap[getPreMappedId(child->getTypeId())].push_back(const_cast<FeatureNode*>(child));
        }
        collectExons(child, exonIdExonMap);
    }
}

/* get the exon mapping version for all exon records.  This is the
 * maximum exon mapping version, or 1 if all exons are passed through.
 */
static int getExonMappingVersion(FeatureNodeVector& exonFeatures)  {
    int maxMappingVersion = 1;
    for (int i = 0; i < exonFeatures.size(); i++) {
        int mappingVersion = getMappingVersion(exonFeatures[i]->getAttrValue(GxfFeature::EXON_ID_ATTR));
        maxMappingVersion = max(mappingVersion, maxMappingVersion);
    }
    return maxMappingVersion;
}

/* Added mapping version numbers to a exon feature */
static void setExonMappingVersion(FeatureNode* exonFeature,
                                  int mappingVersion) {
    const AttrVal* idAttr = exonFeature->getAttr(GxfFeature::EXON_ID_ATTR);
    setMappingVersionInId(exonFeature, idAttr, mappingVersion);
}

/* Added mapping version numbers to a set of exons feeatures with the same id */
static void setExonMappingVersion(const string& exonId,
                                  FeatureNodeVector& exonFeatures,
                                  FeatureNodeVector* prevExonFeatures) {
    int mappingVersion = (prevExonFeatures != NULL)
        ? getExonMappingVersion(*prevExonFeatures) : 1;
    for (int i = 0; i < exonFeatures.size(); i++) {
        if (isRemapped(exonFeatures[i])) {
            setExonMappingVersion(exonFeatures[i], mappingVersion);
        }
    }
}

/* Added mapping version numbers to all exons */
static void setExonsMappingVersions(const FeatureNode* prevGene,
                                    FeatureNode* gene) {
    // exon id scope is gene
    FeatureTreePolish::ExonIdExonMap prevExonIdExonMap, exonIdExonMap;
    if (prevGene != NULL) {
        collectExons(prevGene, prevExonIdExonMap);
    }
    collectExons(gene, exonIdExonMap);
    for (FeatureTreePolish::ExonIdExonMapIter exonIdIter = exonIdExonMap.begin(); exonIdIter != exonIdExonMap.end(); exonIdIter++) {
        FeatureTreePolish::ExonIdExonMapIter prevExonIdIter = prevExonIdExonMap.find(exonIdIter->first);
        FeatureNodeVector* prevExonFeatures = (prevExonIdIter != prevExonIdExonMap.end()) ? &prevExonIdIter->second : NULL;
        setExonMappingVersion(exonIdIter->first, exonIdIter->second, prevExonFeatures);
    }
}

/* Added mapping version numbers to transcripts of a gene */
bool FeatureTreePolish::setTranscriptsMappingVersions(FeatureNode* gene) const {
    bool transcriptsSame = true;
    for (int i = 0; i < gene->getNumChildren(); i++) {
        if (isRemapped(gene->getChild(i))) {
            if (not setTranscriptMappingVersion(gene->getChild(i))) {
                transcriptsSame = false;
            }
        }
    }
    return transcriptsSame;
}


/* Added mapping version numbers */
void FeatureTreePolish::setGeneMappingVersion(FeatureNode* gene) const {
    bool transcriptsSame = setTranscriptsMappingVersions(gene) ;

    const FeatureNode* prevGene = getPrevMappedFeature(gene);
    bool geneSame = transcriptsSame and
        ((prevGene == NULL) or compareGeneFeatures(prevGene, gene));
    int mappingVersion = getFeatureMappingVersion(prevGene, geneSame);
    recursiveSetMappingVersion(gene, GxfFeature::GENE_ID_ATTR, GxfFeature::GENE_HAVANA_ATTR, mappingVersion);
    setExonsMappingVersions(prevGene, gene);
}

/* last minute fix-ups */
void FeatureTreePolish::polishGene(FeatureNode* gene) const {
    // n.b. must renumber exons first, otherwise different exon numbers
    // with previous version will cause false mapping version increments.
    renumberGeneExons(gene);
    if (isRemapped(gene)) {
        setGeneMappingVersion(gene);
    }
}


