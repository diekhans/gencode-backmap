#include "annotationSet.hh"
#include <stdexcept>
#include <iostream>
#include "transMap.hh"

/* add a feature to the location map */
void AnnotationSet::addLocationMap(FeatureNode* featureNode) {
    struct LocationLink* locationLink =  static_cast<struct LocationLink*>(needMem(sizeof(struct LocationLink)));  // zeros memory
    locationLink->featureNode = featureNode;
    genomeRangeTreeAddValList(fLocationMap, toCharStr(featureNode->fFeature->fSeqid),
                              featureNode->fFeature->fStart, featureNode->fFeature->fEnd, locationLink);
}

/* build location map on first use */
void AnnotationSet::buildLocationMap() {
    assert(fLocationMap == NULL);
    fLocationMap = genomeRangeTreeNew();
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        addLocationMap(fGenes[iGene]);
        for (size_t iTrans = 0; iTrans < fGenes[iGene]->fChildren.size(); iTrans++) {
            addLocationMap(fGenes[iGene]->fChildren[iTrans]);
        }
    }
}

/* free the location map tree and data */
void AnnotationSet::freeLocationMap() {
    struct hashCookie chromCookie = hashFirst(fLocationMap->jkhash);
    struct hashEl *chromEl;
    while ((chromEl = hashNext(&chromCookie)) != NULL) {
        struct range *ranges = genomeRangeTreeList(fLocationMap, chromEl->name);
        for (struct range *r = ranges; r != NULL; r = r->next) {
            struct LocationLink *tll, *tlls = static_cast<struct LocationLink*>(r->val);
            while ((tll = static_cast<struct LocationLink*>(slPopHead(&tlls))) != NULL) {
                free(tll);
            }
        }
    }
    genomeRangeTreeFree(&fLocationMap);
}

/* link a gene or transcript feature into the maps */
void AnnotationSet::addFeature(FeatureNode* featureNode) {
    assert(featureNode->isGeneOrTranscript());
    GxfFeature* feature = featureNode->fFeature;
    // record by id and name
    fIdFeatureMap[getBaseId(feature->getTypeId())].push_back(featureNode);
    if (featureNode->getHavanaTypeId() != "") {
        fIdFeatureMap[getBaseId(featureNode->getHavanaTypeId())].push_back(featureNode);
    }
    // save gene/transcript name, although not on small non-coding, as they are
    // not unique.
    if ((feature->getTypeName() != "") && (not featureNode->isAutomaticSmallNonCodingGene())) {
        fNameFeatureMap[feature->getTypeName()].push_back(featureNode);
    }
    if (fLocationMap != NULL) {
        addLocationMap(featureNode);
    }
}

/* add a gene the maps */
void AnnotationSet::addGene(FeatureNode* geneTree) {
    assert(geneTree->isGene());
    fGenes.push_back(geneTree);
    addFeature(geneTree);
    for (size_t i = 0; i < geneTree->fChildren.size(); i++) {
        FeatureNode* transcriptTree = geneTree->fChildren[i];
        if (transcriptTree->fFeature->fType != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcriptTree->fFeature->toString());
        }
        addFeature(transcriptTree);
    }
}

/* process a record, loading into table or discarding */
void AnnotationSet::processRecord(GxfParser *gxfParser,
                                      GxfRecord* gxfRecord) {
    if (instanceOf(gxfRecord, GxfFeature)) {
        GxfFeature* geneFeature = dynamic_cast<GxfFeature*>(gxfRecord);
        FeatureNode* geneTree = GeneTree::geneTreeFactory(gxfParser, geneFeature);
        addGene(geneTree);
    } else {
        delete gxfRecord;
    }
}


/* get a target gene or transcript node from an index by name or id */
FeatureNode* AnnotationSet::getFeatureNodeByKey(const string& key,
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
FeatureNode* AnnotationSet::getFeatureNodeById(const string& id,
                                                   const string& seqIdForParCheck) const {
    return getFeatureNodeByKey(getBaseId(id), fIdFeatureMap, seqIdForParCheck);
}

/* get a target gene or transcript node with same name or NULL.
 * special handling for PARs. Getting node is used if you need whole tree. */
FeatureNode* AnnotationSet::getFeatureNodeByName(const string& name,
                                                 const string& seqIdForParCheck) const {
    return getFeatureNodeByKey(name, fNameFeatureMap, seqIdForParCheck);
}

/* get a target gene or transcript with same base or NULL.
 * special handling for PARs/ */
GxfFeature* AnnotationSet::getFeatureById(const string& id,
                                              const string& seqIdForParCheck) const {
    FeatureNode* featureNode = getFeatureNodeById(id, seqIdForParCheck);
    if (featureNode == NULL) {
        return NULL;
    } else {
        return featureNode->fFeature;
    }
}

/* find overlapping features */
FeatureNodeVector AnnotationSet::findOverlappingFeatures(const string& seqid,
                                                         int start,
                                                         int end) {
    if (fLocationMap == NULL) {
        buildLocationMap();
    }
    FeatureNodeVector overlapping;
    struct range *over, *overs = genomeRangeTreeAllOverlapping(fLocationMap, toCharStr(seqid), start, end);
    for (over = overs; over != NULL; over = over->next) {
        for (struct LocationLink* tll = static_cast<struct LocationLink*>(over->val); tll != NULL; tll = tll->next) {
            overlapping.push_back(tll->featureNode);
        }
    }
    return overlapping;
}

/* is a feature an overlapping gene passing the specified criteria */
bool AnnotationSet::isOverlappingGene(const FeatureNode* geneTree,
                                      const FeatureNode* overlappingFeature,
                                      float minSimilarity,
                                      bool manualOnlyTranscripts) {
    return overlappingFeature->isGene()
        and (geneTree->getMaxTranscriptSimilarity(overlappingFeature,
                                                  manualOnlyTranscripts) >= minSimilarity);
}

/* find overlapping genes with minimum similarity at the transcript level */
FeatureNodeVector AnnotationSet::findOverlappingGenes(const FeatureNode* geneTree,
                                                      float minSimilarity,
                                                      bool manualOnlyTranscripts) {
    FeatureNodeVector overlappingFeatures 
        = findOverlappingFeatures(geneTree->fFeature->fSeqid,
                                  geneTree->fFeature->fStart,
                                  geneTree->fFeature->fEnd);
    FeatureNodeVector overlappingGenes;
    for (int iFeat = 0; iFeat < overlappingFeatures.size(); iFeat++) {
        if (isOverlappingGene(geneTree, overlappingFeatures[iFeat], minSimilarity, manualOnlyTranscripts)) {
            overlappingGenes.push_back(overlappingFeatures[iFeat]);
        }
    }
    return overlappingGenes;
}

/* constructor, load gene and transcript objects from a GxF */
AnnotationSet::AnnotationSet(const string& gxfFile,
                             const GenomeSizeMap* genomeSizes):
    fLocationMap(NULL),
    fGenomeSizes(genomeSizes) {
    GxfParser* gxfParser = GxfParser::factory(gxfFile);
    GxfRecord* gxfRecord;
    while ((gxfRecord = gxfParser->next()) != NULL) {
        processRecord(gxfParser, gxfRecord);
    }
    delete gxfParser;
}

/* destructor */
AnnotationSet::~AnnotationSet() {
    if (fLocationMap != NULL) {
        freeLocationMap();
    }
    for (int i = 0; i <fGenes.size(); i++) {
        delete fGenes[i];
    }
}


/* write a sequence region record */
void AnnotationSet::outputSeqRegion(const string& seqId,
                                    int size,
                                    GxfWriter& gxfFh) {
    gxfFh.write(string("##sequence-region ") + seqId + " 1 " + toString(size));
}

/* output GFF3 mapped ##sequence-region if not already written */
void AnnotationSet::outputMappedSeqRegionIfNeed(const FeatureNode* geneTree,
                                                GxfWriter& gxfFh) {
    if (gxfFh.getFormat() == GFF3_FORMAT) {
        const string& seqId = geneTree->fFeature->fSeqid;
        if (fGenomeSizes->have(seqId) and (not checkRecordSeqRegionWritten(seqId))) {
            outputSeqRegion(seqId, fGenomeSizes->get(seqId), gxfFh);
        }
    }
}

/*
 * recursive output of a GxF feature tree
 */
void AnnotationSet::outputFeature(const FeatureNode* featureNode,
                                  GxfWriter& gxfFh) const {
    gxfFh.write(featureNode->fFeature);
    for (size_t i = 0; i < featureNode->fChildren.size(); i++) {
        outputFeature(featureNode->fChildren[i], gxfFh);
    }
}

/* output genes */
void AnnotationSet::write(GxfWriter& gxfFh) {
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        outputMappedSeqRegionIfNeed(fGenes[iGene], gxfFh);
        outputFeature(fGenes[iGene], gxfFh);
    }
}


