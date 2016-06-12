/**
 * Handling mapping of a transcript
 */
#ifndef transcriptMapper_hh
#define transcriptMapper_hh
class TransMap;
class PslMapping;
class FeatureTransMap;
class Feature;
class AnnotationSet;
#include "feature.hh"
#include "resultFeatures.hh"
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
    const Feature* fTargetGene;                     // target annotations for this transcript, if any, to help
    const Feature* fTargetTranscript;               // selecting between multiple mappings.
    static const bool debug = 0;
    
    static FeatureVector getExons(const Feature* transcript);
    PslMapping* allExonsTransMap(const Feature* transcript) const;
    static const TransMapVector makeViaExonsTransMap(const PslMapping* exonsMapping);
    PslMapping* featurePslMap(const Feature* feature);
    TransMappedFeature mapFeature(const Feature* feature);
    TransMappedFeature mapFeatures(const Feature* feature);
    ResultFeatures mapTranscriptFeature(const Feature* transcript);

    public:
    /* constructor, targetAnnotations can be NULL */
    TranscriptMapper(const TransMap* genomeTransMap,
                     const Feature* transcript,
                     const AnnotationSet* targetAnnotations,
                     bool srcSeqInMapping,
                     ostream* transcriptPslFh);

    /* destructor */
    ~TranscriptMapper();
    
    /*
     * map one transcript's annotations.  Fill in transcript
     */
    ResultFeatures mapTranscriptFeatures(const Feature* transcript);
};

#endif
