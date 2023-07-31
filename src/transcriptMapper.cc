#include "transcriptMapper.hh"
#include "featureMapper.hh"
#include "featureTree.hh"
#include "pslOps.hh"
#include "pslMapping.hh"
#include "frame.hh"
#include "featureTransMap.hh"
#include "typeOps.hh"
#include "annotationSet.hh"
#include <stdexcept>
#include <iostream>

static const bool debug = 0;


/* build transcript exons PSL to query and mapping to target genome.
 * Return NULL if no mappings for whatever reason.*/
PslMapping* TranscriptMapper::allExonsTransMap(const FeatureNode* transcript) const {
    const string& qName(transcript->getAttr(GxfFeature::TRANSCRIPT_ID_ATTR)->getVal());
    FeatureNodeVector exons;
    transcript->getMatchingType(exons, GxfFeature::EXON);

    // get alignment of exons to srcGenome and to targetGenome
    PslMapping* exonsMapping = FeatureTransMap(fGenomeTransMap).mapFeatures(qName, exons);
    if (exonsMapping == NULL) {
        return NULL; // not in mapping chains
    }
        
    // filter for same chrom and sort using more evidence
    exonsMapping->filterSameTarget();
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
    transMaps.push_back(TransMap::factoryFromPsl(exonsMapping->getSrcPsl(), true)); // swap map genomeA to exons
    transMaps.push_back(TransMap::factoryFromPsl(exonsMapping->getBestMappedPsl(), false)); // exons to genomeB
    return transMaps;
}

/* get PSL of feature mapping */
PslMapping* TranscriptMapper::featurePslMap(const FeatureNode* feature) {
    const AttrVal* idAttr = feature->findAttr(GxfFeature::ID_ATTR);
    const string& featureId = (idAttr != NULL) ? idAttr->getVal() : "someFeature";
    PslMapping* pslMapping = fViaExonsFeatureTransMap->mapFeature(featureId, feature);
    pslMapping->sortMappedPsls();
    return pslMapping;
}

/* map one feature, linking in a child feature. */
TransMappedFeature TranscriptMapper::mapFeature(const FeatureNode* feature) {
    PslMapping* pslMapping = (fViaExonsFeatureTransMap != NULL) ? featurePslMap(feature) : NULL;
    TransMappedFeature transMappedFeature = FeatureMapper::map(feature, pslMapping);
    transMappedFeature.setRemapStatus(fSrcSeqInMapping);
    delete pslMapping;
    return transMappedFeature;
}

/* recursive map features below transcript */
TransMappedFeature TranscriptMapper::mapFeatures(const FeatureNode* feature) {
    TransMappedFeature transMappedFeature = mapFeature(feature);
    for (int iChild = 0; iChild < feature->getNumChildren(); iChild++) {
        TransMappedFeature childFeatures = mapFeatures(feature->getChild(iChild));
        FeatureMapper::updateParents(transMappedFeature, childFeatures);
    }
    return transMappedFeature;
}

/* create a new transcript record that covers the alignment */
ResultFeatureTrees TranscriptMapper::mapTranscriptFeature(const FeatureNode* transcript) {
    ResultFeatureTrees mappedTranscript(transcript);
    struct psl* mappedPsl = (fExonsMapping != NULL) ? fExonsMapping->getBestMappedPsl() : NULL;
    if (mappedPsl != NULL) {
        // transcript for mapped PSLs
        mappedTranscript.mapped
            = FeatureMapper::mapBounding(transcript,
                                         string(mappedPsl->tName),
                                         mappedPsl->tStart, mappedPsl->tEnd,
                                         charToString(pslTStrand(mappedPsl)));
        mappedTranscript.mapped->setNumMappings(fExonsMapping->getMappedPsls().size());
    }

    // if any parts was unmapped, also need a copy of the original transcript
    if ((mappedPsl == NULL) or (not pslQueryFullyMapped(mappedPsl))) {
        mappedTranscript.unmapped = FeatureMapper::mapBounding(transcript);
    }
    return mappedTranscript;
}

/* constructor, targetAnnotations can be NULL */
TranscriptMapper::TranscriptMapper(const TransMap* genomeTransMap,
                                   const FeatureNode* transcript,
                                   const AnnotationSet* targetAnnotations,
                                   bool srcSeqInMapping,
                                   ostream* transcriptPslFh):
    fGenomeTransMap(genomeTransMap), 
    fSrcSeqInMapping(srcSeqInMapping),
    fExonsMapping(NULL),
    fViaExonsFeatureTransMap(NULL),
    fTargetGene(NULL),
    fTargetTranscript(NULL) {
    assert(transcript->getType() == GxfFeature::TRANSCRIPT);

    // if available, find target transcripts to use in selecting multiple mappings.  PAR is handled
    // my only allowing same chromosome.
    if (targetAnnotations != NULL) {
        fTargetGene = targetAnnotations->getFeatureById(transcript->getAttrValue(GxfFeature::GENE_ID_ATTR),
                                                        transcript->fFeature->getSeqid());
        fTargetTranscript = targetAnnotations->getFeatureById(transcript->getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR),
                                                              transcript->fFeature->getSeqid());
    }

    // map all exons together, this will be used to project the other features
    fExonsMapping = allExonsTransMap(transcript);
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
 * map one transcript's annotations.  Fill in transcript
 */
ResultFeatureTrees TranscriptMapper::mapTranscriptFeatures(const FeatureNode* transcript) {
    // project features via exons (including redoing exons)
    ResultFeatureTrees mappedTranscript = mapTranscriptFeature(transcript);
    TransMappedFeature mappedTranscriptSet(mappedTranscript);
    for (int iChild = 0; iChild < transcript->getNumChildren(); iChild++) {
        TransMappedFeature transMappedFeature = mapFeatures(transcript->getChild(iChild));
        FeatureMapper::updateParents(mappedTranscriptSet, transMappedFeature);
    }

    mappedTranscript.setBoundingFeatureRemapStatus(fSrcSeqInMapping);
    mappedTranscript.setNumMappingsAttr();
    return mappedTranscript;
}
