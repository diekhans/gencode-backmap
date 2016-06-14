#include "annotationSet.hh"
#include <stdexcept>
#include <iostream>
#include "transMap.hh"
#include "featureIO.hh"
#include "gxfIO.hh"

// FIXME: should writer in featureIO


/* add a feature to the location map */
void AnnotationSet::addLocationMap(Feature* feature) {
    struct LocationLink* locationLink =  static_cast<struct LocationLink*>(needMem(sizeof(struct LocationLink)));  // zeros memory
    locationLink->feature = feature;
    genomeRangeTreeAddValList(fLocationMap, toCharStr(feature->getSeqid()),
                              feature->getStart(), feature->getEnd(), locationLink);
}

/* build location map on first use */
void AnnotationSet::buildLocationMap() {
    assert(fLocationMap == NULL);
    fLocationMap = genomeRangeTreeNew();
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        addLocationMap(fGenes[iGene]);
        for (size_t iTrans = 0; iTrans < fGenes[iGene]->getChildren().size(); iTrans++) {
            addLocationMap(fGenes[iGene]->getChild(iTrans));
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
void AnnotationSet::addFeature(Feature* feature) {
    assert(feature->isGeneOrTranscript());
    // record by id and name
    fIdFeatureMap[getBaseId(feature->getTypeId())].push_back(feature);
    if (feature->getHavanaTypeId() != "") {
        fIdFeatureMap[getBaseId(feature->getHavanaTypeId())].push_back(feature);
    }
    // save gene/transcript name, although not on small non-coding, as they are
    // not unique.
    if ((feature->getTypeName() != "") && (not feature->isAutomaticSmallNonCodingGene())) {
        fNameFeatureMap[feature->getTypeName()].push_back(feature);
    }
    if (fLocationMap != NULL) {
        addLocationMap(feature);
    }
}

/* add a gene the maps */
void AnnotationSet::addGene(Feature* gene) {
    assert(gene->isGene());
    fGenes.push_back(gene);
    addFeature(gene);
    for (size_t i = 0; i < gene->getChildren().size(); i++) {
        Feature* transcript = gene->getChild(i);
        if (transcript->getType() != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcript->toString());
        }
        addFeature(transcript);
    }
}

/* get a target gene or transcript node from an index by name or id */
Feature* AnnotationSet::getFeatureByKey(const string& key,
                                                const FeatureMap& featureMap,
                                                const string& seqIdForParCheck) const {
    FeatureMapConstIter it = featureMap.find(key);
    if (it == featureMap.end()) {
        return NULL;
    } else if (it->second.size() == 2) {
        if (it->second[0]->getSeqid() == seqIdForParCheck) {
            return it->second[0];
        } else if (it->second[1]->getSeqid() == seqIdForParCheck) {
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
Feature* AnnotationSet::getFeatureById(const string& id,
                                               const string& seqIdForParCheck) const {
    return getFeatureByKey(getBaseId(id), fIdFeatureMap, seqIdForParCheck);
}

/* get a target gene or transcript node with same name or NULL.
 * special handling for PARs. Getting node is used if you need whole tree. */
Feature* AnnotationSet::getFeatureByName(const string& name,
                                                 const string& seqIdForParCheck) const {
    return getFeatureByKey(name, fNameFeatureMap, seqIdForParCheck);
}

/* find overlapping features */
FeatureVector AnnotationSet::findOverlappingFeatures(const string& seqid,
                                                     int start,
                                                     int end) {
    if (fLocationMap == NULL) {
        buildLocationMap();
    }
    FeatureVector overlapping;
    struct range *over, *overs = genomeRangeTreeAllOverlapping(fLocationMap, toCharStr(seqid), start, end);
    for (over = overs; over != NULL; over = over->next) {
        for (struct LocationLink* tll = static_cast<struct LocationLink*>(over->val); tll != NULL; tll = tll->next) {
            overlapping.push_back(tll->feature);
        }
    }
    return overlapping;
}

/* is a feature an overlapping gene passing the specified criteria */
bool AnnotationSet::isOverlappingGene(const Feature* gene,
                                      const Feature* overlappingFeature,
                                      float minSimilarity,
                                      bool manualOnlyTranscripts) {
    return overlappingFeature->isGene()
        and (gene->getMaxTranscriptSimilarity(overlappingFeature,
                                                  manualOnlyTranscripts) >= minSimilarity);
}

/* find overlapping genes with minimum similarity at the transcript level */
FeatureVector AnnotationSet::findOverlappingGenes(const Feature* gene,
                                                      float minSimilarity,
                                                      bool manualOnlyTranscripts) {
    FeatureVector overlappingFeatures 
        = findOverlappingFeatures(gene->getSeqid(),
                                  gene->getStart(),
                                  gene->getEnd());
    FeatureVector overlappingGenes;
    for (int iFeat = 0; iFeat < overlappingFeatures.size(); iFeat++) {
        if (isOverlappingGene(gene, overlappingFeatures[iFeat], minSimilarity, manualOnlyTranscripts)) {
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
    FeatureParser parser(gxfFile);
    Feature* gene;
    while ((gene = parser.nextGene()) != NULL) {
        addGene(gene);
    }
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
void AnnotationSet::outputMappedSeqRegionIfNeed(const Feature* gene,
                                                GxfWriter& gxfFh) {
    if (gxfFh.getFormat() == GFF3_FORMAT) {
        const string& seqId = gene->getSeqid();
        if (fGenomeSizes->have(seqId) and (not checkRecordSeqRegionWritten(seqId))) {
            outputSeqRegion(seqId, fGenomeSizes->get(seqId), gxfFh);
        }
    }
}

/*
 * recursive output of a GxF feature tree
 */
void AnnotationSet::outputFeature(const Feature* feature,
                                  GxfWriter& gxfFh) const {
    gxfFh.write(feature);
    for (size_t i = 0; i < feature->getChildren().size(); i++) {
        outputFeature(feature->getChild(i), gxfFh);
    }
}

/* print for debugging */
void AnnotationSet::dump(ostream& fh) const {
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        fGenes[iGene]->dump(fh);
    }
}

/* output genes */
void AnnotationSet::write(GxfWriter& gxfFh) {
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        outputMappedSeqRegionIfNeed(fGenes[iGene], gxfFh);
        outputFeature(fGenes[iGene], gxfFh);
    }
}
