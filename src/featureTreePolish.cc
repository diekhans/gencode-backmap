/*
 * Post mapping corrections to feature tree.
 */
#include "featureTreePolish.hh"
#include <stdexcept>
#include <algorithm>
#include "annotationSet.hh"
#include <iostream>
#include "frame.hh"
#include "globals.hh"

static const bool DEBUG = false;

/* Gap closing in mappings:
 * Looked at fixing in mapping alignments, however there can be some complex alignment
 * indel patterns that result in relatively simple gapped mappings, so we fix after the
 * fact instead.
 */
static const int maxCloseGapSize = 6;

/* node, feature and index in transcript */
class NodeIndex {
    public:
    FeatureNode* feat;
    int idx;

    NodeIndex(FeatureNode* feature,
              int idx):
        feat(feature), idx(idx) {
    }

    NodeIndex():
        feat(NULL), idx(-1) {
    }

    bool isNull() const {
        return feat == NULL;
    }
};

typedef vector<NodeIndex> NodeIndexVector;


/* sort features in all transcripts */
static void sortGeneContaining(FeatureNode* gene) {
    for (int i = 0; i < gene->getNumChildren(); i++) {
        gene->getChild(i)->getChildren().sortContaining();
    }
}

    
/* replace exon_number attr */
static void setExonNumber(FeatureNode *feature, int exonNum) {
    feature->getAttrs().update(AttrVal(GxfFeature::EXON_NUMBER_ATTR, toString(exonNum)));
}

/* next exon feature index, starting at position */
static NodeIndex getNextExon(FeatureNode* transcript,
                             int iFeat) {
    // assumes exons are in order
    for (int jFeat = iFeat; jFeat < transcript->getNumChildren(); jFeat++) {
        if (transcript->getChild(jFeat)->getType() == GxfFeature::EXON) {
            return NodeIndex(transcript->getChild(jFeat), jFeat);
        }
    }
    return NodeIndex();
}

/* Make sure features is after the exon in the transcript array, as this
 * ordering is counted on by where the merged features get inserted. */
static void ensureFeatureOrder(const NodeIndex& exon,
                               int iFeat) {
    if (iFeat < exon.idx) {
        throw logic_error("ensureFeatureOrder: feature order is not as expected");
    }
}

/* CDS feature index overlapping exon  */
static NodeIndex getOverlappingCds(FeatureNode* transcript,
                                   const NodeIndex& exon) {
    // assumes CDS are after exon, but doesn't insist on immediately after
    for (int iFeat = exon.idx + 1; iFeat < transcript->getNumChildren(); iFeat++) {
        if ((transcript->getChild(iFeat)->getType() == GxfFeature::CDS)
            and (transcript->getChild(iFeat)->overlaps(exon.feat))) {
            ensureFeatureOrder(exon, iFeat);
            return NodeIndex(transcript->getChild(iFeat), iFeat);
        }
    }
    return NodeIndex();
}

/*  features overlapping exon, including CDS. */
static NodeIndexVector getExonOverlapping(FeatureNode* transcript,
                                          const NodeIndex& exon) {
    // assumes overlap are after exon, but doesn't insist on immediately after
    NodeIndexVector feats;
    for (int iFeat = exon.idx + 1; iFeat < transcript->getNumChildren(); iFeat++) {
        if ((transcript->getChild(iFeat)->getType() != GxfFeature::EXON)
            and (transcript->getChild(iFeat)->overlaps(exon.feat))) {
            ensureFeatureOrder(exon, iFeat);
            feats.push_back(NodeIndex(transcript->getChild(iFeat), iFeat));
        }
    }
    return feats;
}

/* Merge two feature records */
static string getMergePhase(FeatureNode* feature1, FeatureNode* feature2) {
   if (feature1->getPhase() != ".") {
        if (feature1->getStrand() == "+") {
            return feature1->getPhase();
        } else {
            return feature2->getPhase();
        }
   } else {
       return ".";
   }
}

/* FeatureNode with new feature and existing attributes */
static FeatureNode* mkMergeNode(FeatureNode* transcript,
                                GxfFeature* newFeature, FeatureNode* feature) {
    return new FeatureNode(newFeature, transcript, feature->fRemapStatus, feature->fTargetStatus, feature->fNumMappings);
}

/* Merge two feature records */
static FeatureNode* mergeFeatureRecs(FeatureNode* transcript,
                                     FeatureNode* feature1, FeatureNode* feature2) {
    if (gVerbose) {
        cerr << "merge: " << transcript->getTypeId() << ": " << feature1->getType()
             << ": " << feature1->getStart() << "-" << feature1->getEnd()
             << " with " << feature2->getStart() << "-" << feature2->getEnd()
             << endl;
    }
    if (feature1->getType() != feature2->getType()) {
        throw logic_error("mergeFeatureRecs: feature types not the same");
    }
    if ((feature1->getNumChildren() > 0) or (feature2->getNumChildren() > 0)) {
        throw logic_error("mergeFeatureRecs: assumption of no transcript grandchildren is not true");
    }
    GxfFeature* newFeature = new GxfFeature(feature1->getSeqid(), feature1->getSource(), feature1->getType(),
                                            feature1->getStart(), feature2->getEnd(), feature1->getScore(),
                                            feature1->getStrand(), getMergePhase(feature1, feature2), feature1->getAttrs());
    return mkMergeNode(transcript, newFeature, feature1);
}

/* find feature of specified type or -1 */
static int findOverType(const string& featType,
                        NodeIndexVector& overs,
                        int startIdx = 0) {
    for (int i = 0; i < overs.size(); i++) {
        if (overs[i].feat->getType() == featType) {
            return i;
        }
    }
    return -1;
}

/* merge if possible */
static FeatureNode* mergeCopyOther(FeatureNode* transcript,
                                   FeatureNode* feature,
                                   NodeIndexVector& overs2,
                                   BoolVector& done2) {
    int i2 = findOverType(feature->getType(), overs2);
    if (i2 < 0) {
        return feature->cloneFeature();
    } else {
        FeatureNode* newFeat = mergeFeatureRecs(transcript, feature, overs2[i2].feat);
        done2[i2] = true;
        return newFeat;
    }
}

/* merge or copy overlapping other features */
static FeatureNodeVector mergeCopyOthers(FeatureNode* transcript,
                                         NodeIndexVector& overs1,
                                         NodeIndexVector& overs2) {
    FeatureNodeVector merged;
    BoolVector done2;
    done2.resize(overs2.size());

    for (int i1 = 0; i1 < overs1.size(); i1++) {
        FeatureNode* feat = mergeCopyOther(transcript, overs1[i1].feat, overs2, done2);
        merged.push_back(feat);
    }

    // copy those that haven't been merged
    for (int i2 = 0; i2 < overs2.size(); i2++) {
        if (not done2[i2]) {
            FeatureNode* feat = overs2[i2].feat->cloneFeature();
            merged.push_back(feat);
        }
    }
    return merged;
}

/* clear out nodes that have been replaced */
static void removeReplacedNodes(FeatureNode* transcript,
                                NodeIndexVector& overs) {
    for (int i = 0; i < overs.size(); i++) {
        delete transcript->removeChild(transcript->getChildIdx(overs[i].feat));
    }
}

/* insert or append a node */
static void insertNode(FeatureNode* transcript,
                       int iFeat,
                       FeatureNode* feat) {
    FeatureNodeVector& children = transcript->getChildren();
    assert(iFeat <= children.size());
    if (iFeat >= children.size()) {
        children.push_back(feat);
    } else {
        children.insert(children.begin() + (iFeat + 1), feat);
    }
}
                             
/* add new nodes at position iFeat */
static void addNewNodes(FeatureNode* transcript,
                        int iFeat,
                        FeatureNode* mergedExon,
                        FeatureNodeVector& mergedOthers) {
    insertNode(transcript, iFeat++, mergedExon);
    for (int i = 0; i < mergedOthers.size(); i++) {
        insertNode(transcript, iFeat++, mergedOthers[i]);
    }
}

/* find the insert point for replaced records after deletions */
static int findInsertPoint(const FeatureNode* transcript,
                           const FeatureNode* mergedExon) {
    const FeatureNodeVector& children = transcript->getChildren();
    for (int iChild = 0; iChild < children.size(); iChild++) {
        if (children[iChild]->getStart0() >= mergedExon->getEnd()) {
            return iChild;
        }
    }
    return children.size();
}

/* Merge two exon records, index in feature array of exon1 stays the same,
 * everything above it might change */
static void mergeExonRecs(FeatureNode* transcript,
                          NodeIndex& exon1, NodeIndex& exon2) {
    // N.B. this recreates records even if they could just moved to keep book
    // keeping easier, as this is not performance-sensitive.
    FeatureNode* mergedExon = mergeFeatureRecs(transcript, exon1.feat, exon2.feat);

    // other features
    NodeIndexVector overs1 = getExonOverlapping(transcript, exon1);
    NodeIndexVector overs2 = getExonOverlapping(transcript, exon2);
    FeatureNodeVector mergedOthers = mergeCopyOthers(transcript, overs1, overs2);

    // delete replaced features
    delete transcript->removeChild(transcript->getChildIdx(exon1.feat));
    delete transcript->removeChild(transcript->getChildIdx(exon2.feat));
    removeReplacedNodes(transcript, overs1);
    removeReplacedNodes(transcript, overs2);
    
    addNewNodes(transcript, findInsertPoint(transcript, mergedExon),
                mergedExon, mergedOthers);
}

static bool preservesFrame(NodeIndex& cds1,
                           NodeIndex& cds2) {
    int gapSize = cds2.feat->getStart0() - cds1.feat->getEnd();
    Frame frame1(Frame::fromPhaseStr(cds1.feat->getPhase()));
    Frame frame2(Frame::fromPhaseStr(cds2.feat->getPhase()));
    if (cds1.feat->getStrand() == "+") {
        return frame1.incr(cds1.feat->length() + gapSize) == frame2;
    } else {
        return frame2.incr(cds2.feat->length() + gapSize) == frame1;
    }
}

/* Check and possible merge two exon */
static bool canMergeExonRecs(FeatureNode* transcript,
                             NodeIndex& exon1,
                             NodeIndex& exon2) {
    if ((exon2.feat->getStart0() - exon1.feat->getEnd()) > maxCloseGapSize) {
        return false;   // gap over threshold
    }
    NodeIndex cds1 = getOverlappingCds(transcript, exon1);
    NodeIndex cds2 = getOverlappingCds(transcript, exon2);
    if ((not cds1.isNull()) and (not cds2.isNull())
        and (not preservesFrame(cds1, cds2))) {
        return false;   // CDS on both sides and not keeping frame
    }
    return true;
}

/* Merge gaps in exons causes by gaps in alignments. */
static void mergeTranscriptGaps(FeatureNode* transcript) {
    int iFeat = 0;
    bool anyMerged = false;
    while (true) {
        NodeIndex exon1 = getNextExon(transcript, iFeat);
        if (exon1.isNull()) {
            break;
        }
        iFeat = exon1.idx + 1;
        NodeIndex exon2 = getNextExon(transcript, iFeat);
        if (exon2.isNull()) {
            break;
        }
        if (canMergeExonRecs(transcript, exon1, exon2)) {
            mergeExonRecs(transcript, exon1, exon2);
            anyMerged = true;
        } else {
            iFeat++;
        }
    }
    if (anyMerged) {
        transcript->getChildren().sortContaining();
    }
}

/* Merge gaps in all . */
static void mergeGeneGaps(FeatureNode* gene) {
    for (int i = 0; i < gene->getNumChildren(); i++) {
        mergeTranscriptGaps(gene->getChild(i));
    }
}

/* count of exon features */
static int countNumExons(FeatureNode* transcript) {
    int cnt = 0;
    for (int i = 0; i < transcript->getNumChildren(); i++) {
        if (transcript->getChild(i)->getType() == GxfFeature::EXON) {
            cnt++;
        }
    }
    return cnt;
}

/* renumber all exons in a transcript */
static void renumberTranscriptExons(FeatureNode* transcript) {
    // logic here assumes exons come before features contained
    // in them
    if (transcript->getStrand() == "+") {
        int exonNum = 0;
        for (int i = 0; i < transcript->getNumChildren(); i++) {
            if (transcript->getChild(i)->getType() == GxfFeature::EXON) {
                exonNum++;
            }
            setExonNumber(transcript->getChild(i), exonNum);
        }
    } else {
        int exonNum = countNumExons(transcript);
        for (int i = 0; i < transcript->getNumChildren(); i++) {
            setExonNumber(transcript->getChild(i), exonNum);
            if (transcript->getChild(i)->getType() == GxfFeature::EXON) {
                exonNum--;
            }
        }
    }
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
        prevChildren.sortContaining();
        FeatureNodeVector newChildren(newParent->getChildren());
        newChildren.sortContaining();

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
    sortGeneContaining(gene);
    if (isRemapped(gene)) {
        // n.b. must renumber exons first, otherwise different exon numbers
        // with previous version will cause false mapping version increments.
        mergeGeneGaps(gene);
        renumberGeneExons(gene);
        setGeneMappingVersion(gene);
    }
}


