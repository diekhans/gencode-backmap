#include "transcriptMapper.hh"
#include "featureMapper.hh"
#include "featureTree.hh"
#include "jkinclude.hh"
#include "pslOps.hh"
#include "gxf.hh"
#include "pslMapping.hh"
#include "frame.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include "targetAnnotations.hh"
#include <stdexcept>
#include <iostream>

/* get exon features */
GxfFeatureVector TranscriptMapper::getExons(const FeatureNode* transcriptTree) {
    GxfFeatureVector exons;
    transcriptTree->getMatching(exons, [](const GxfFeature* f) {
            return f->fType == GxfFeature::EXON;
        });
    return exons;
}

/* build transcript exons PSL to query and mapping to target genome.
 * Return NULL if no mappings for whatever reason.*/
PslMapping* TranscriptMapper::allExonsTransMap(const FeatureNode* transcriptTree) const {
    const string& qName(transcriptTree->fFeature->getAttr(GxfFeature::TRANSCRIPT_ID_ATTR)->fVal);
    GxfFeatureVector exons = getExons(transcriptTree);
    // get alignment of exons to srcGenome and to targetGenome
    PslMapping* exonsMapping = FeatureTransMap(fGenomeTransMap).mapFeatures(qName, exons);
    if (exonsMapping == NULL) {
        return NULL;  // source sequence not in map; FIXME: not sure this should even get here in this case
    }
    // resort using more evidence
    exonsMapping->sortMappedPsls(fTargetTranscript, fTargetGene);
    if (debug) {
        exonsMapping->dump(cerr, "Transcript Exons:", "    ");
    }
    if (not exonsMapping->haveMappings()) {
        delete exonsMapping;
        return NULL;
    } else {
        return exonsMapping;
    }
}

/* create transMap objects used to do two level mapping via exons. */
const TransMapVector TranscriptMapper::makeViaExonsTransMap(const PslMapping* exonsMapping) {
    TransMapVector transMaps;
    transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fSrcPsl, true)); // swap map genomeA to exons
    transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fMappedPsl, false)); // exons to genomeB
    return transMaps;
}

/* create a new transcript record that covers the alignment */
void TranscriptMapper::mapTranscriptFeature(FeatureNode* transcriptNode) {
    struct psl* mappedPsl = fExonsMapping->fMappedPsl;
    // transcript for mapped PSLs
    FeatureMapper::mapBounding(transcriptNode,
                               string(mappedPsl->tName),
                               mappedPsl->tStart, mappedPsl->tEnd,
                               charToString(pslTStrand(mappedPsl)));
    transcriptNode->fNumMappings = fExonsMapping->fMappedPsls.size();
    // if any parts was unmapped, also need a copy of the original transcript
    if (not pslQueryFullyMapped(mappedPsl)) {
        FeatureMapper::mapBounding(transcriptNode);
    }
}

/* recursive map features below transcript */
void TranscriptMapper::mapFeatures(FeatureNode* featureNode) {
    const AttrVal* idAttr = featureNode->fFeature->findAttr(GxfFeature::ID_ATTR);
    const string& nodeId = (idAttr != NULL) ? idAttr->fVal : "someFeature";
    PslMapping* pslMapping = fViaExonsFeatureTransMap->mapFeature(nodeId, featureNode->fFeature);
    FeatureMapper::map(featureNode, pslMapping);
    for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
       mapFeatures(featureNode->fChildren[iChild]);
    }
    delete pslMapping;
}

/* do work of mapping features when we have transmap mapping alignments
 * and know something will map */
void TranscriptMapper::mapTranscriptFeaturesViaExons(FeatureNode* transcriptTree) {
    mapTranscriptFeature(transcriptTree);
    for (int iChild = 0; iChild < transcriptTree->fChildren.size(); iChild++) {
        mapFeatures(transcriptTree->fChildren[iChild]);
    }
}

/* recursive process features that are unmapped */
void TranscriptMapper::processUnmappedFeatures(FeatureNode* featureNode) {
    FeatureMapper::map(featureNode, NULL);
    for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
        processUnmappedFeatures(featureNode->fChildren[iChild]);
    }
}

/* constructor, targetAnnotations can be NULL */
TranscriptMapper::TranscriptMapper(const TransMap* genomeTransMap,
                                   FeatureNode* transcriptTree,
                                   const TargetAnnotations* targetAnnotations,
                                   bool srcSeqInMapping,
                                   ostream* transcriptPslFh):
    fGenomeTransMap(genomeTransMap), 
    fSrcSeqInMapping(srcSeqInMapping),
    fExonsMapping(NULL),
    fViaExonsFeatureTransMap(NULL),
    fTargetGene(NULL),
    fTargetTranscript(NULL) {
    assert(transcriptTree->fFeature->fType == GxfFeature::TRANSCRIPT);

    // if available, find target transcripts to use in selecting multiple mappings.  Special handling
    // for PAR requires sequence id.
    if (targetAnnotations != NULL) {
        fTargetGene = targetAnnotations->getFeature(transcriptTree->fFeature->getAttrValue(GxfFeature::GENE_ID_ATTR),
                                                    transcriptTree->fFeature->fSeqid);
        fTargetTranscript = targetAnnotations->getFeature(transcriptTree->fFeature->getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR),
                                                          transcriptTree->fFeature->fSeqid);
    }

    // map all exons together, this will be used to project the other exons
    fExonsMapping = allExonsTransMap(transcriptTree);
    if (fExonsMapping != NULL) {
        if (transcriptPslFh != NULL) {
            fExonsMapping->writeMapped(*transcriptPslFh);
        }
        fVaiExonsTransMaps = makeViaExonsTransMap(fExonsMapping);
        fViaExonsFeatureTransMap = new FeatureTransMap(fVaiExonsTransMaps);
    }
}

/* destructor */
TranscriptMapper::~TranscriptMapper() {
    delete fExonsMapping;
    fVaiExonsTransMaps.free();
    delete fViaExonsFeatureTransMap;
}

/*
 * map one transcript's annotations.  Fill in transcriptTree
 */
void TranscriptMapper::mapTranscriptFeatures(FeatureNode* transcriptTree) {
    // project features via exons (including redoing exons)
    if (fViaExonsFeatureTransMap != NULL) {
        mapTranscriptFeaturesViaExons(transcriptTree);
    } else {
        processUnmappedFeatures(transcriptTree);
    }
    transcriptTree->recursiveCalcRemapStatus(fSrcSeqInMapping);
    FeatureMapper::updateIds(transcriptTree);
    transcriptTree->setNumMappingsAttr();
}
