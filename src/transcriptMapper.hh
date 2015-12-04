/**
 * Handling mapping of a transcript
 */
#ifndef transcriptMapper_hh
#define transcriptMapper_hh
class TransMap;
class PslMapping;
class FeatureTransMap;
class GxfFeature;
class TargetAnnotations;
#include "featureTree.hh"
#include "gxf.hh"
#include "transMap.hh"

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
    
    static GxfFeatureVector getExons(const FeatureNode* transcriptTree);
    PslMapping* allExonsTransMap(const FeatureNode* transcriptTree) const;
    static const TransMapVector makeViaExonsTransMap(const PslMapping* exonsMapping);
    PslMapping* featurePslMap(const FeatureNode* featureNode);
    TransMappedFeature mapNodeFeature(const FeatureNode* featureNode);
    TransMappedFeature mapFeatures(const FeatureNode* featureNode);
    ResultFeatureTrees mapTranscriptFeature(const FeatureNode* transcriptNode);

    public:
    /* constructor, targetAnnotations can be NULL */
    TranscriptMapper(const TransMap* genomeTransMap,
                     const FeatureNode* transcriptTree,
                     const TargetAnnotations* targetAnnotations,
                     bool srcSeqInMapping,
                     ostream* transcriptPslFh);

    /* destructor */
    ~TranscriptMapper();
    
    /*
     * map one transcript's annotations.  Fill in transcriptTree
     */
    ResultFeatureTrees mapTranscriptFeatures(const FeatureNode* transcriptTree);
};

#endif
