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
    const string& qName(transcriptTree->fFeature->getAttr(GxfFeature::TRANSCRIPT_ID_ATTR)->getVal());
    GxfFeatureVector exons = getExons(transcriptTree);
    // get alignment of exons to srcGenome and to targetGenome
    PslMapping* exonsMapping = FeatureTransMap(fGenomeTransMap).mapFeatures(qName, exons);
    if (exonsMapping == NULL) {
        return NULL;  // source sequence not in map
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

/* get PSL of feature mapping */
PslMapping* TranscriptMapper::featurePslMap(const FeatureNode* featureNode) {
    const AttrVal* idAttr = featureNode->fFeature->findAttr(GxfFeature::ID_ATTR);
    const string& nodeId = (idAttr != NULL) ? idAttr->getVal() : "someFeature";
    return fViaExonsFeatureTransMap->mapFeature(nodeId, featureNode->fFeature);
}

/* map one nodes feature, linking in a child nodes. */
TransMappedFeature TranscriptMapper::mapNodeFeature(const FeatureNode* featureNode) {
    PslMapping* pslMapping = (fViaExonsFeatureTransMap != NULL) ? featurePslMap(featureNode) : NULL;
    TransMappedFeature transMappedFeature = FeatureMapper::map(featureNode, pslMapping);
    transMappedFeature.setRemapStatus(fSrcSeqInMapping);
    delete pslMapping;
    return transMappedFeature;
}

/* recursive map features below transcript */
TransMappedFeature TranscriptMapper::mapFeatures(const FeatureNode* featureNode) {
    TransMappedFeature transMappedFeature = mapNodeFeature(featureNode);
    for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
        TransMappedFeature childFeatureNodes = mapFeatures(featureNode->fChildren[iChild]);
        FeatureMapper::updateParents(transMappedFeature, childFeatureNodes);
    }
    return transMappedFeature;
}

/* create a new transcript record that covers the alignment */
ResultFeatureTrees TranscriptMapper::mapTranscriptFeature(const FeatureNode* transcriptNode) {
    ResultFeatureTrees mappedTranscript(transcriptNode);
    struct psl* mappedPsl = (fExonsMapping != NULL) ? fExonsMapping->fMappedPsl : NULL;
    if (mappedPsl != NULL) {
        // transcript for mapped PSLs
        mappedTranscript.mapped
            = FeatureMapper::mapBounding(transcriptNode,
                                         string(mappedPsl->tName),
                                         mappedPsl->tStart, mappedPsl->tEnd,
                                         charToString(pslTStrand(mappedPsl)));
        mappedTranscript.mapped->fNumMappings = fExonsMapping->fMappedPsls.size();
    }

    // if any parts was unmapped, also need a copy of the original transcript
    if ((mappedPsl == NULL) or (not pslQueryFullyMapped(mappedPsl))) {
        mappedTranscript.unmapped = FeatureMapper::mapBounding(transcriptNode);
    }
    return mappedTranscript;
}

/* constructor, targetAnnotations can be NULL */
TranscriptMapper::TranscriptMapper(const TransMap* genomeTransMap,
                                   const FeatureNode* transcriptTree,
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
        fTargetGene = targetAnnotations->getFeatureById(transcriptTree->fFeature->getAttrValue(GxfFeature::GENE_ID_ATTR),
                                                        transcriptTree->fFeature->fSeqid);
        fTargetTranscript = targetAnnotations->getFeatureById(transcriptTree->fFeature->getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR),
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
ResultFeatureTrees TranscriptMapper::mapTranscriptFeatures(const FeatureNode* transcriptTree) {
    // project features via exons (including redoing exons)
    ResultFeatureTrees mappedTranscript = mapTranscriptFeature(transcriptTree);
    TransMappedFeature mappedTranscriptSet(mappedTranscript);
    for (int iChild = 0; iChild < transcriptTree->fChildren.size(); iChild++) {
        TransMappedFeature transMappedFeature = mapFeatures(transcriptTree->fChildren[iChild]);
        FeatureMapper::updateParents(mappedTranscriptSet, transMappedFeature);
    }

    mappedTranscript.setBoundingFeatureRemapStatus(fSrcSeqInMapping);
    mappedTranscript.setNumMappingsAttr();
    return mappedTranscript;
}
