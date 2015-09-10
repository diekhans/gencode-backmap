/**
 * Handling mapping of a transcript
 */
#ifndef transcriptMapper_hh
#define transcriptMapper_hh
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


/**
 * Class to map a single transcript and subfeatures
 * All mapping is done by a two-level alignment.
 * with mapping alignments.
 *  - to map cdsA in  annotation of genomeA to genomeB:
 *    - an alignment of [genomeA->exonsA]
 *    - an alignment of [exonsA->genomeB]
 *    - do a two level transmap:
 *      cdsA->genomeA =>  [genomeA->exonsA] => cdsA->exonA => [exonsA->genomeB] => cdsA->genomeB
 */
class TranscriptMapper {
    private:
    const TransMap* fGenomeTransMap;
    const bool fSrcSeqInMapping;                 // do we have source sequence in genomic mapps
    const PslMapping* fExonsMapping;            // exons as psl and genome mapping of exons.
    TransMapVector fVaiExonsTransMaps;          // transmap objects that are combined
    const FeatureTransMap* fViaExonsFeatureTransMap;   // two-level transmap, NULL if can't map (owned)
    const GxfFeature* fTargetGene;                     // target annotations for this transcript, if any, to help
    const GxfFeature* fTargetTranscript;               // selecting between multiple mappings.
    static const bool debug = 0;
    
    /* get exon features */
    static GxfFeatureVector getExons(const FeatureNode* transcriptTree) {
        GxfFeatureVector exons;
        transcriptTree->getMatching(exons, [](const GxfFeature* f) {
                return f->fType == GxfFeature::EXON;
            });
        return exons;
    }

    /* build transcript exons PSL to query and mapping to target genome.
     * Return NULL if no mappings for whatever reason.*/
    PslMapping* allExonsTransMap(const FeatureNode* transcriptTree) const {
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
    static const TransMapVector makeViaExonsTransMap(const PslMapping* exonsMapping) {
        TransMapVector transMaps;
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fSrcPsl, true)); // swap map genomeA to exons
        transMaps.push_back(TransMap::factoryFromPsls(exonsMapping->fMappedPsl, false)); // exons to genomeB
        return transMaps;
    }

    /* create a new transcript record that covers the alignment */
    void mapTranscriptFeature(FeatureNode* transcriptNode) {
        struct psl* mappedPsl = fExonsMapping->fMappedPsl;
        // transcript for mapped PSLs
        FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping,
                                   string(mappedPsl->tName),
                                   mappedPsl->tStart, mappedPsl->tEnd,
                                   charToString(pslQStrand(mappedPsl)));
        transcriptNode->fNumMappings = fExonsMapping->fMappedPsls.size();
        // if any parts was unmapped, also need a copy of the original transcript
        if (not pslQueryFullyMapped(mappedPsl)) {
            FeatureMapper::mapBounding(transcriptNode, fSrcSeqInMapping);
        }
    }
    
    /* recursive map features below transcript */
    void mapFeatures(FeatureNode* featureNode) {
        const AttrVal* idAttr = featureNode->fFeature->findAttr(GxfFeature::ID_ATTR);
        const string& nodeId = (idAttr != NULL) ? idAttr->fVal : "someFeature";
        PslMapping* pslMapping = fViaExonsFeatureTransMap->mapFeature(nodeId, featureNode->fFeature);
        FeatureMapper::map(featureNode, pslMapping, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
           mapFeatures(featureNode->fChildren[iChild]);
        }
        delete pslMapping;
    }
    
    /* do work of mapping features when we have transmap mapping alignments
     * and know something will map */
    void mapTranscriptFeaturesViaExons(FeatureNode* transcriptTree) {
        mapTranscriptFeature(transcriptTree);
        for (int iChild = 0; iChild < transcriptTree->fChildren.size(); iChild++) {
            mapFeatures(transcriptTree->fChildren[iChild]);
        }
    }

    /* recursive process features that are unmapped */
    void processUnmappedFeatures(FeatureNode* featureNode) {
        FeatureMapper::map(featureNode, NULL, fSrcSeqInMapping);
        for (int iChild = 0; iChild < featureNode->fChildren.size(); iChild++) {
            processUnmappedFeatures(featureNode->fChildren[iChild]);
        }
    }

    public:
    /* constructor, targetAnnotations can be NULL */
    TranscriptMapper(const TransMap* genomeTransMap,
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

        // if available, find target transcripts to use in selecting multiple mappings
        if (targetAnnotations != NULL) {
            fTargetGene = targetAnnotations->get(transcriptTree->fFeature->getAttrValue(GxfFeature::GENE_ID_ATTR));
            fTargetTranscript = targetAnnotations->get(transcriptTree->fFeature->getAttrValue(GxfFeature::TRANSCRIPT_ID_ATTR));
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
    ~TranscriptMapper() {
        delete fExonsMapping;
        fVaiExonsTransMaps.free();
        delete fViaExonsFeatureTransMap;
    }
    
    /*
     * map one transcript's annotations.  Fill in transcriptTree
     */
    void mapTranscriptFeatures(FeatureNode* transcriptTree) {
        // project features via exons (including redoing exons)
        if (fViaExonsFeatureTransMap != NULL) {
            mapTranscriptFeaturesViaExons(transcriptTree);
        } else {
            processUnmappedFeatures(transcriptTree);
        }
        FeatureMapper::updateIds(transcriptTree);
        transcriptTree->setNumMappingsAttr();
    }
};

#endif
