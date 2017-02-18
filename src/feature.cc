/* 
 * GENCODE features as a tree.
 */
#include "feature.hh"
#include <algorithm>
#include <iostream>
#include "typeOps.hh"

/* Remap status attribute name */
const string REMAP_STATUS_ATTR = "remap_status";

/* Attribute name used for original id before remap */
const string REMAP_ORIGINAL_ID_ATTR = "remap_original_id";

/* Attribute name used for original localization before remap */
const string REMAP_ORIGINAL_LOCATION_ATTR = "remap_original_location";

/* Attribute name for count of mappings, set on transcripts or genes */
const string REMAP_NUM_MAPPINGS_ATTR = "remap_num_mappings";

/* Attribute name for target of mapping */
const string REMAP_TARGET_STATUS_ATTR = "remap_target_status";

/* Attribute indicating target gene was substituted due to   */
const string REMAP_SUBSTITUTED_MISSING_TARGET_ATTR = "remap_substituted_missing_target";

/* ensembl non-coding gene biotypes to skip */
static const char* automaticNonCodingGeneBiotypes[] = {
    "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "ribozyme", "rRNA", "scaRNA",
    "snoRNA", "snRNA", "sRNA", NULL
};

/* compare chrom names to emulate GENCODE sorting */
static bool chromLessThan(const string& a, const string& b) {
    // chrom vs non-chrom; ucsc names have chr_accession, so check for that too
    bool aIsChr = (a.find("chr") == 0) or (a.find("_") < 0);
    bool bIsChr = (b.find("chr") == 0) or (b.find("_") < 0);
    if (aIsChr and (not bIsChr)) {
        return true;
    } else if (bIsChr and (not aIsChr)) {
        return false;
    } else if ((not aIsChr) and (not bIsChr)) {
        return aIsChr < bIsChr; // not a chrom
    }
    // autosomes, or X,Y, or M
    bool aIsAuto = isdigit(a[3]);
    bool bIsAuto = isdigit(b[3]);
    if (aIsAuto and (not bIsAuto)) {
        return true;
    } else if (bIsAuto and (not aIsAuto)) {
        return false;
    }
    
    // put chrM last
    bool aIsChrM = (a == "chrM");
    bool bIsChrM = (b == "chrM");
    if (aIsChrM and (not bIsChrM)) {
        return false;
    } else if (bIsChrM and (not aIsChrM)) {
        return true;
    } else if (aIsChrM and bIsChrM) {
        return false;
    }
    // both chroms
    if (a.size() != b.size()) {
        return a.size() < b.size();  // chr10 vs chr1
    }
    return a < b;
}

/* sort the vector in a predictable order.  This is not necessary what
 * will be in the GxF file by GENCODE conventions. */
void FeatureVector::sort() {
    std::sort(begin(), end(),
              [](const Feature* a, const Feature* b) -> bool {
                  if (a->getSeqid() != b->getSeqid()) {
                      return chromLessThan(a->getSeqid(), b->getSeqid());
                  } else if (a->getStart() != b->getStart()) {
                      return a->getStart() < b->getStart();
                  } else {
                      return a->getEnd() < b->getEnd();
                  }
              });
}


/* is ensembl small non-coding gene */
bool Feature::isAutomaticSmallNonCodingGene() const {
    if (fSource != "ENSEMBL") {
        return false;
    }
    const string& bioType = getTypeBiotype();
    for (int i = 0; automaticNonCodingGeneBiotypes[i] != NULL; i++) {
        if (bioType == automaticNonCodingGeneBiotypes[i]) {
            return true;
        }
    }
    return false;
}

/* set the remap number of mappings attribute on this node.  not recursive,
 * since it's only set on gene/transcript */
void Feature::setNumMappingsAttr() {
    getAttrs().update(AttrVal(REMAP_NUM_MAPPINGS_ATTR, ::toString(fNumMappings)));
}

/* recursively set the remap status attribute */
void Feature::rsetRemapStatusAttr() {
    if (fRemapStatus != REMAP_STATUS_NONE) {
        getAttrs().update(AttrVal(REMAP_STATUS_ATTR, remapStatusToStr(fRemapStatus)));
    }
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetRemapStatusAttr();
    }
}

/* recursively set the target status attribute */
void Feature::rsetTargetStatusAttr() {
    if (isGeneOrTranscript()) {
        if (fTargetStatus != TARGET_STATUS_NA) {
            getAttrs().update(AttrVal(REMAP_TARGET_STATUS_ATTR, targetStatusToStr(fTargetStatus)));
        }
        for (int i = 0; i < fChildren.size(); i++) {
            fChildren[i]->rsetTargetStatusAttr();
        }
    }
}

/* recursively set the target status attribute on fFeature node. */
void Feature::rsetSubstitutedMissingTargetAttr(const string& targetVersion) {
    getAttrs().update(AttrVal(REMAP_SUBSTITUTED_MISSING_TARGET_ATTR, targetVersion));
    for (int i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetSubstitutedMissingTargetAttr(targetVersion);
    }
}

/* depth-first output */
void Feature::write(ostream& fh) const {
    fh << toString() << endl;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->write(fh);
    }
}


/* recursively set the remap status. */
void Feature::rsetRemapStatus(RemapStatus remapStatus) {
    fRemapStatus = remapStatus;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetRemapStatus(remapStatus);
    }
}

/* Set the target status. Not recursive */
void Feature::setTargetStatus(TargetStatus targetStatus) {
    fTargetStatus = targetStatus;
}

/* Recursively set the target status. */
void Feature::rsetTargetStatus(TargetStatus targetStatus) {
    fTargetStatus = targetStatus;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->rsetTargetStatus(targetStatus);
    }
}

/* do any child belond to the specified status */
bool Feature::anyChildWithRemapStatus(unsigned remapStatusSet) const {
    for (size_t i = 0; i < fChildren.size(); i++) {
        if ((fChildren[i]->fRemapStatus & remapStatusSet) != 0) {
            return true;
        }
    }
    return false;
}

/* do all child have belong to the specified status set */
bool Feature::allChildWithRemapStatus(unsigned remapStatusSet) const {
    for (size_t i = 0; i < fChildren.size(); i++) {
        if ((fChildren[i]->fRemapStatus & remapStatusSet) == 0) {
            return false;
        }
    }
    return true;
}

/* get the size of a transcript, in exons */
int Feature::getTranscriptExonSize() const {
    assert(isTranscript());
    int size = 0;
    for (int iExon = 0; iExon < fChildren.size(); iExon++) {
        if (fChildren[iExon]->isExon()) {
            size += fChildren[iExon]->size();
        }
    }
    return size;
}

/* count overlapping bases */
int Feature::getOverlapAmount(const Feature* other) const {
    int maxStart = max(other->fStart, fStart);
    int minEnd = min(other->fEnd, fEnd);
    return (maxStart <= minEnd) ? (minEnd-maxStart)+1 : 0;
}

/* count exon overlap with exons of another transcript */
int Feature::countExonOverlap(const Feature* exon1,
                                  const Feature* trans2) const {
    int totalOverlap = 0;
    for (int iFeat2 = 0; iFeat2 < trans2->fChildren.size(); iFeat2++) {
        if (trans2->fChildren[iFeat2]->isExon()) {
            const Feature* exon2 = trans2->fChildren[iFeat2];
            totalOverlap += exon1->getOverlapAmount(exon2);
        }
    }
    return totalOverlap;
}

/* get exon similarity */
float Feature::getExonSimilarity(const Feature* trans2) const {
    assert(isTranscript());
    assert(trans2->isTranscript());
    int totalOverlap = 0;
    for (int iFeat1 = 0; iFeat1 < fChildren.size(); iFeat1++) {
        if (fChildren[iFeat1]->isExon()) {
            totalOverlap += countExonOverlap(fChildren[iFeat1], trans2);
        }
    }
    return float(2*totalOverlap)/float(getTranscriptExonSize() + trans2->getTranscriptExonSize());
}

/* get the maximum transcript similarity for a gene and a transcript  */
float Feature::getMaxTranscriptSimilarity(const Feature* gene2,
                                              const Feature* trans1,
                                              bool manualOnlyTranscripts) const {
    float maxSimilarity = 0.0;
    for (int iTrans2 = 0; (iTrans2 < gene2->fChildren.size()) && (maxSimilarity < 1.0); iTrans2++) {
        const Feature* trans2 = gene2->fChildren[iTrans2];
        if ((not manualOnlyTranscripts) or (not trans2->isAutomatic())) {
            float similarity = trans1->getExonSimilarity(trans2);
            maxSimilarity = max(similarity, maxSimilarity);
       }
    }
    return maxSimilarity;
}

/* get the maximum transcript similarity for a gene */
float Feature::getMaxTranscriptSimilarity(const Feature* gene2,
                                              bool manualOnlyTranscripts) const {
    assert(isGene());
    assert(gene2->isGene());
    float maxSimilarity = 0.0;
    for (int iTrans1 = 0; (iTrans1 < fChildren.size()) && (maxSimilarity < 1.0); iTrans1++) {
        const Feature* trans1 = fChildren[iTrans1];
        if ((not manualOnlyTranscripts) or (not trans1->isAutomatic())) {
            maxSimilarity = max(maxSimilarity,
                                getMaxTranscriptSimilarity(gene2, trans1, manualOnlyTranscripts));
        }
    }
    return maxSimilarity;
}
/* clone tree */
Feature* Feature::cloneTree() const {
    Feature *newFeature = static_cast<Feature*>(clone());
    for (int i = 0; i < fChildren.size(); i++) {
        newFeature->fChildren.push_back(static_cast<Feature*>(fChildren[i]->cloneTree()));
    }
    return newFeature;
}

/* print node for debugging */
/* recursively print for debugging */
void Feature::dump(ostream& fh) const {
    fh << toString() << endl;
    for (size_t i = 0; i < fChildren.size(); i++) {
        fChildren[i]->dump(fh);
    }
}

/* function to create a new feature based on type */
Feature* featureFactory(const string& seqid, const string& source, const string& type,
                        int start, int end, const string& score, const string& strand,
                        const string& phase, const AttrVals& attrs) {
    return new Feature(seqid, source, type, start, end, score, strand, phase, attrs);
}
