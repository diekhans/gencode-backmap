#include "targetAnnotations.hh"
#include <stdexcept>
#include <iostream>

/* link a gene or transcript feature into the maps */
void TargetAnnotations::loadFeature(FeatureNode* featureNode) {
    assert(featureNode->isGeneOrTranscript());
    GxfFeature* feature = featureNode->fFeature;
    // record by id and name
    string baseId = getBaseId(feature->getTypeId());
    fIdFeatureMap[baseId].push_back(featureNode);
    if (feature->getTypeName() != "") {
        fNameFeatureMap[feature->getTypeName()].push_back(featureNode);
    }
    
    struct TargetLocationLink* locationLink =  static_cast<struct TargetLocationLink*>(needMem(sizeof(struct TargetLocationLink)));  // zeros memory
    locationLink->featureNode = featureNode;
    genomeRangeTreeAddValList(fLocationMap, toCharStr(feature->fSeqid),
                              feature->fStart, feature->fEnd, locationLink);
}

/* link a gene or transcript feature into the maps */
void TargetAnnotations::loadGene(FeatureNode* geneTree) {
    fGenes.push_back(geneTree);
    loadFeature(geneTree);
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        loadFeature(transcriptTree);
    }
}

/* process a record, loading into table or discarding */
void TargetAnnotations::processRecord(GxfParser *gxfParser,
                                      GxfRecord* gxfRecord) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        GxfFeature* geneFeature = dynamic_cast<GxfFeature*>(gxfRecord);
        FeatureNode* geneTree = GeneTree::geneTreeFactory(gxfParser, geneFeature);
        loadGene(geneTree);
    } else {
        delete gxfRecord;
    }
}


/* get a target gene or transcript node from an index by name or id */
FeatureNode* TargetAnnotations::getFeatureNodeByKey(const string& key,
                                                    const FeatureMap& featureMap,
                                                    const string& seqIdForParCheck) const {
    FeatureMapConstIter it = featureMap.find(key);
    if (it == featureMap.end()) {
        return NULL;
    } else if (it->second.size() == 2) {
        if (it->second[0]->fFeature->fSeqid == seqIdForParCheck) {
            return it->second[0];
        } else if (it->second[1]->fFeature->fSeqid == seqIdForParCheck) {
            return it->second[1];
        } else {
            throw logic_error("PAR target feature hack confused: " + key);
        }
    } else if (it->second.size() > 2) {
        throw logic_error("too many nodes for key: " + key);
    } else {
        return it->second[0];
    }
}

/* get a target gene or transcript node with same base id or NULL.
 * special handling for PARs. Getting node is used if you need whole tree. */
FeatureNode* TargetAnnotations::getFeatureNodeById(const string& id,
                                                   const string& seqIdForParCheck) const {
    return getFeatureNodeByKey(getBaseId(id), fIdFeatureMap, seqIdForParCheck);
}

/* get a target gene or transcript node with same name or NULL.
 * special handling for PARs. Getting node is used if you need whole tree. */
FeatureNode* TargetAnnotations::getFeatureNodeByName(const string& name,
                                                     const string& seqIdForParCheck) const {
    return getFeatureNodeByKey(name, fNameFeatureMap, seqIdForParCheck);
}

/* get a target gene or transcript with same base or NULL.
 * special handling for PARs/ */
GxfFeature* TargetAnnotations::getFeatureById(const string& id,
                                              const string& seqIdForParCheck) const {
    FeatureNode* featureNode = getFeatureNodeById(id, seqIdForParCheck);
    if (featureNode == NULL) {
        return NULL;
    } else {
        return featureNode->fFeature;
    }
}

/* find overlapping features */
FeatureNodeVector TargetAnnotations::findOverlappingFeatures(const string& seqid,
                                                             int start,
                                                             int end) {
    FeatureNodeVector overlapping;
    struct range *over, *overs = genomeRangeTreeAllOverlapping(fLocationMap, toCharStr(seqid), start, end);
    for (over = overs; over != NULL; over = over->next) {
        for (struct TargetLocationLink* tll = static_cast<struct TargetLocationLink*>(over->val); tll != NULL; tll = tll->next) {
            overlapping.push_back(tll->featureNode);
        }
    }
    return overlapping;
}

/* constructor, load gene and transcript objects from a GxF */
TargetAnnotations::TargetAnnotations(const string& gxfFile):
    fLocationMap(genomeRangeTreeNew()) {
    GxfParser* gxfParser = GxfParser::factory(gxfFile);
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        processRecord(gxfParser, gxfRecord);
    }
    delete gxfParser;
}

/* free the location map tree and data */
void TargetAnnotations::freeLocationMap() {
    struct hashCookie chromCookie = hashFirst(fLocationMap->jkhash);
    struct hashEl *chromEl;
    while ((chromEl = hashNext(&chromCookie)) != NULL) {
        struct range *ranges = genomeRangeTreeList(fLocationMap, chromEl->name);
        for (struct range *r = ranges; r != NULL; r = r->next) {
            struct TargetLocationLink *tll, *tlls = static_cast<struct TargetLocationLink*>(r->val);
            while ((tll = static_cast<struct TargetLocationLink*>(slPopHead(&tlls))) != NULL) {
                free(tll);
            }
        }
    }
    genomeRangeTreeFree(&fLocationMap);
}

/* destructor */
TargetAnnotations::~TargetAnnotations() {
    freeLocationMap();
    for (int i = 0; i <fGenes.size(); i++) {
        delete fGenes[i];
    }
}

