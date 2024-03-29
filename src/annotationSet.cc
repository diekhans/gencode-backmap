#include "annotationSet.hh"
#include <stdexcept>
#include <iostream>
#include "transMap.hh"
#include "globals.hh"
#include "gxf.hh"

// FIXME: should writer in featureIO


/* add a feature to the location map */
void AnnotationSet::addLocationMap(FeatureNode* feature) {
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
        for (size_t iTrans = 0; iTrans < fGenes[iGene]->getNumChildren(); iTrans++) {
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

/* generate key with PAR */
string AnnotationSet::mkFeatureIdKey(const string& typeId,
                                     const string& chrom) const {
    return chrom + "/" + typeId;
}

/* insert feature in a map */
void AnnotationSet::insertInFeatureMap(const string& key,
                                       FeatureNode* feature,
                                       FeatureMap& featureMap) {
    if (featureMap.find(key) != featureMap.end()) {
        if (gVerbose) {
            cerr << "NOTE: key already in FeatureMap: " << key << endl;
        }
    }
    featureMap[key].push_back(feature);
}

/* check if we already have this feature with a different version */
void AnnotationSet::checkForMultiIdVersions(FeatureNode* feature) const {
    const FeatureNodeVector& features = getFeaturesById(getBaseId(feature->getTypeId()));
    for (int i = 0; i < features.size(); i++) {
        if (features[i]->fFeature->getTypeId() != feature->getTypeId()) {
            cerr << "Note: feature id " + feature->getTypeId() + " already exists with different version: " + features[i]->fFeature->getTypeId() << endl;
        }
    }
}


/* link a gene or transcript feature into the maps */
void AnnotationSet::addFeature(FeatureNode* feature) {
    if (not fWarnIdDiffVersions) {
        checkForMultiIdVersions(feature);
    }
    
    assert(feature->isGeneOrTranscript());
    // record by base id and name with chrom
    insertInFeatureMap(mkFeatureIdKey(getBaseId(feature->getTypeId()), feature->fFeature->getSeqid()),
                       feature, fIdFeatureChromMap);
    if (feature->getHavanaTypeId() != "") {
        insertInFeatureMap(mkFeatureIdKey(getBaseId(feature->getHavanaTypeId()), feature->fFeature->getSeqid()),
                           feature, fIdFeatureChromMap);
    }
    // save gene name when real and unique
    if (feature->isGene() and useGeneNameForMappingKey(feature)) {
        insertInFeatureMap(mkFeatureIdKey(feature->getTypeName(), feature->fFeature->getSeqid()),
                           feature, fNameFeatureChromMap);
    }
    if (fLocationMap != NULL) {
        addLocationMap(feature);
    }

    // record for checking for ids of different versions
    insertInFeatureMap(getBaseId(feature->getTypeId()), feature, fIdFeatureMap);
}

/* add a gene the maps */
void AnnotationSet::addGene(FeatureNode* gene) {
    assert(gene->isGene());
    fGenes.push_back(gene);
    addFeature(gene);
    for (size_t i = 0; i < gene->getNumChildren(); i++) {
        FeatureNode* transcript = gene->getChild(i);
        if (transcript->getType() != GxfFeature::TRANSCRIPT) {
            throw logic_error("gene record has child that is not of type transcript: " + transcript->toString());
        }
        addFeature(transcript);
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

/* get a target gene or transcript node from an index by name or id.
 * if the name or id is duplicated in the GxF, it can't be used as an index and
 * NULL is returned.
 */
FeatureNode* AnnotationSet::getFeatureByKey(const string& baseId,
                                            const string& chrom,
                                            const FeatureMap& featureMap) const {
    string key = mkFeatureIdKey(baseId, chrom);
    FeatureMapConstIter it = featureMap.find(key);
    if (it == featureMap.end()) {
        return NULL;
    } else if (it->second.size() > 1) {
        return NULL;
    } else {
        return it->second[0];
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
            overlapping.push_back(tll->feature);
        }
    }
    return overlapping;
}

/* is a feature an overlapping gene passing the specified criteria */
bool AnnotationSet::isOverlappingGene(const FeatureNode* gene,
                                      const FeatureNode* overlappingFeature,
                                      float minSimilarity,
                                      bool manualOnlyTranscripts) {
    return overlappingFeature->isGene()
        and (gene->getMaxTranscriptSimilarity(overlappingFeature,
                                                  manualOnlyTranscripts) >= minSimilarity);
}

/* find overlapping genes with minimum similarity at the transcript level */
FeatureNodeVector AnnotationSet::findOverlappingGenes(const FeatureNode* gene,
                                                      float minSimilarity,
                                                      bool manualOnlyTranscripts) {
    FeatureNodeVector overlappingFeatures 
        = findOverlappingFeatures(gene->getSeqid(),
                                  gene->getStart(),
                                  gene->getEnd());
    FeatureNodeVector overlappingGenes;
    for (int iFeat = 0; iFeat < overlappingFeatures.size(); iFeat++) {
        if (isOverlappingGene(gene, overlappingFeatures[iFeat], minSimilarity, manualOnlyTranscripts)) {
            overlappingGenes.push_back(overlappingFeatures[iFeat]);
        }
    }
    return overlappingGenes;
}

/* constructor, load gene and transcript objects from a GxF */
AnnotationSet::AnnotationSet(const string& gxfFile,
                             const GenomeSizeMap* genomeSizes,
                             bool warnIdDiffVersions):
    fWarnIdDiffVersions(warnIdDiffVersions),
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
    for (int i = 0; i < fGenes.size(); i++) {
        delete fGenes[i];
    }
}

/* sort gene in gencode order */
static void sortGencodeGene(FeatureNode* gene) {
    gene->getChildren().sortChromAttrName(GxfFeature::TRANSCRIPT_NAME_ATTR); // sorts transcripts
    for (int i = 0; i < gene->getNumChildren(); i++) {
        gene->getChild(i)->getChildren().sortContainingStrand();
    }
}

void AnnotationSet::sortGencode() {
    sortGenes();
    for (int i = 0; i < fGenes.size(); i++) {
        sortGencodeGene(fGenes[i]);
    }
}

/* write a sequence region record */
void AnnotationSet::outputSeqRegion(const string& seqId,
                                    int size,
                                    GxfWriter& gxfFh) {
    gxfFh.write(string("##sequence-region ") + seqId + " 1 " + toString(size));
}

/* output GFF3 mapped ##sequence-region if not already written */
void AnnotationSet::outputMappedSeqRegionIfNeed(const FeatureNode* gene,
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
void AnnotationSet::outputFeature(const FeatureNode* feature,
                                  GxfWriter& gxfFh) const {
    gxfFh.write(feature->getGxfFeature());
    for (size_t i = 0; i < feature->getNumChildren(); i++) {
        outputFeature(feature->getChild(i), gxfFh);
    }
}

/* print genes for debugging */
void AnnotationSet::dump(ostream& fh) const {
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        fGenes[iGene]->dump(fh);
    }
}

/* dump one of the id/name maps */
void AnnotationSet::dumpFeatureMap(const FeatureMap& featureMap,
                                   const string& label,
                                   ostream& fh) const {
    fh << ">>> " << label << endl;
    for (FeatureMapConstIter it = featureMap.begin(); it != featureMap.end(); it++) {
        fh << it->first << ":";
        for (int i = 0; i < it->second.size(); i++) {
            fh << " " << it->second[i]->getTypeId();
        }
        fh << endl;
    }
}

/* print id maps for debugging */
void AnnotationSet::dumpIdMaps(ostream& fh) const {
    dumpFeatureMap(fIdFeatureChromMap, "IdFeatureChromMap", fh);
    dumpFeatureMap(fNameFeatureChromMap, "NameFeatureChromMap", fh);
}

/* output genes */
void AnnotationSet::write(GxfWriter& gxfFh) {
    for (int iGene = 0; iGene < fGenes.size(); iGene++) {
        outputMappedSeqRegionIfNeed(fGenes[iGene], gxfFh);
        outputFeature(fGenes[iGene], gxfFh);
    }
}
