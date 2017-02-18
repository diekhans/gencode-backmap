/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxfRecord.hh"
#include "feature.hh"
#include "typeOps.hh"
#include "resultFeatures.hh"
#include <set>
class TransMap;
class PslMapping;
struct psl;
class PslCursor;
class AnnotationSet;
class BedMap;
class FeatureTreePolish;
class GxfWriter;

/* class that maps a gene to the new assemble */
class GeneMapper {
    public:
    // flags for use target instead of map
    enum {
        useTargetForAutoNonCoding = 0x01,
        useTargetForAutoGenes     = 0x02,
        useTargetForPseudoGenes   = 0x04,
        useTargetForPatchRegions  = 0x08
    };
    private:
    const AnnotationSet* fSrcAnnotations; // source annotations
    const TransMap* fGenomeTransMap;  // genomic mapping
    const AnnotationSet* fTargetAnnotations; // targeted genes/transcripts, maybe NULL
    const AnnotationSet* fPreviousMappedAnotations; // previous version
    const BedMap* fTargetPatchMap; // location of patch regions in target genome
    const string fSubstituteTargetVersion;  // pass through targets when gene new gene doesn't map
    unsigned fUseTargetFlags;  // what targets to force.
    bool fOnlyManualForTargetSubstituteOverlap;  // only check manual transcripts when checking target/map overlap

    /* set of base ids (gene, transcript, havana) and gene names that have been
     * mapped.  Used to prevent output of target genes types that are not being
     * mapped (automatic genes), when type or source changes are already mapped.
     * N.B. Can't use AnnotationSet to track this, due to PAR mappings needing
     * to be mapped twice. */
    StringSet fMappedIdsNames;
    
    int fCurrentGeneNum;  /* used by output info log to logically group features together,
                           * increments each time a gene is process */ 
    
    void outputInfoHeader(ostream& mappingInfoFh) const;
    void outputInfo(const string& recType,
                    const string& featType,
                    const Feature* feature,
                    RemapStatus mappingStatus,
                    int mappingCount,
                    TargetStatus targetStatus,
                    ostream& mappingInfoFh) const;
    void outputSrcGeneInfo(const ResultFeatures* mappedGene,
                           ostream& mappingInfoFh) const;
    void outputMappedGeneInfo(const ResultFeatures* mappedGene,
                              ostream& mappingInfoFh) const;
    void outputUnmappedGeneInfo(const ResultFeatures* mappedGene,
                                ostream& mappingInfoFh) const;
    void outputTargetGeneInfo(const ResultFeatures* mappedGene,
                              const string& targetAction,
                              ostream& mappingInfoFh) const;
    string featureDesc(const Feature* feature) const;
    bool isSrcSeqInMapping(const Feature* feature) const;
    void debugRecordMapped(const Feature* feature,
                           const string& desc,
                           const string& key = "") const;
    void recordGeneMapped(const Feature* gene);
    void recordTranscriptMapped(const Feature* transcript);
    bool checkGeneMapped(const Feature* gene) const ;
    bool checkTranscriptMapped(const Feature* transcript) const;
    bool checkGeneTranscriptsMapped(const Feature* gene) const;
    ResultFeatures processTranscript(const Feature* transcript,
                                     ostream* transcriptPslFh) const;
    ResultFeaturesVector processTranscripts(const Feature* gene,
                                            ostream* transcriptPslFh) const;
    Feature* findMatchingBoundingFeature(const FeatureVector& features,
                                         const Feature* srcFeature) const;
    void copyMappingMetadata(const Feature* origFeature,
                             Feature* newFeature) const;
    void copyGeneMetadata(const Feature* origGene,
                          Feature* newGene) const;
    void forceToUnmapped(ResultFeatures* mappedGene) const;
    void forceToUnmappedDueToRemapStatus(ResultFeatures* mappedGene,
                                         RemapStatus remapStatus) const;
    void forceToUnmappedDueToTargetStatus(ResultFeatures* mappedGene,
                                          TargetStatus targetStatus) const;
    bool hasMixedMappedSeqStrand(const ResultFeatures* mappedGene) const;
    bool hasTargetStatusNonOverlap(const ResultFeatures* mappedGene) const;
    bool hasExcessiveSizeChange(const ResultFeatures* mappedGene) const;

    void setNumGeneMappings(Feature* mappedGeneTree) const;
    bool checkForPathologicalGeneRename(const ResultFeatures* mappedGene,
                                        const Feature* targetGene) const;
    bool shouldSubstituteTarget(const ResultFeatures* mappedGene) const;
    void substituteTarget(ResultFeatures* mappedGene);
    void updateMappedGeneBounds(const Feature* mappedTranscript,
                                string& seqid, string& strand,
                                int& start, int& end) const;
    Feature* buildMappedGeneFeature(const Feature* srcGeneTree,
                                        ResultFeaturesVector& mappedTranscripts) const;
    Feature* buildUnmappedGeneFeature(const Feature* srcGeneTree,
                                          ResultFeaturesVector& mappedTranscripts) const;
    ResultFeatures buildGeneFeature(const Feature* srcGeneTree,
                                        ResultFeaturesVector& mappedTranscripts) const;
    void saveMapped(ResultFeatures& mappedGene,
                    AnnotationSet& mappedSet);
    void saveUnmapped(ResultFeatures& mappedGene,
                      AnnotationSet& unmappedSet);
    const Feature* getTargetAnnotation(const Feature* feature) const;
    TargetStatus getTargetAnnotationStatus(const ResultFeatures* mappedFeature) const;
    const string& getTargetAnnotationBiotype(const ResultFeatures* mappedFeature) const;
    void processGeneLevelMapping(ResultFeatures* mappedGene);
    void setGeneLevelMappingAttributes(ResultFeatures* mappedGene);
    void mapGene(const Feature* srcGeneTree,
                 AnnotationSet& mappedSet,
                 AnnotationSet& unmappedSet,
                 FeatureTreePolish& featureTreePolish,
                 ostream& mappingInfoFh,
                 ostream* transcriptPslFh);
    RemapStatus getNoMapRemapStatus(const Feature* gene) const;
    bool shouldMapGeneType(const Feature* gene) const;
    bool inTargetPatchRegion(const Feature* targetGene);
    bool checkTargetOverlappingMapped(const Feature* targetGene,
                                      AnnotationSet& mappedSet);
    bool shouldIncludeTargetGene(const Feature* gene,
                                 AnnotationSet& mappedSet);
    void copyTargetGene(const Feature* targetGene,
                        AnnotationSet& mappedSet,
                        ostream& mappingInfoFh);
    void copyTargetGenes(AnnotationSet& mappedSet,
                         ostream& mappingInfoFh);
    public:
    /* Constructor */
    GeneMapper(const AnnotationSet* srcAnnotations,
               const TransMap* genomeTransMap,
               const AnnotationSet* targetAnnotations,
               const AnnotationSet* previousMappedAnnotations,
               const BedMap* targetPatchMap,
               const string& substituteTargetVersion,
               unsigned useTargetFlags,
               bool onlyManualForTargetSubstituteOverlap):
        fSrcAnnotations(srcAnnotations),
        fGenomeTransMap(genomeTransMap),
        fTargetAnnotations(targetAnnotations),
        fPreviousMappedAnotations(previousMappedAnnotations),
        fTargetPatchMap(targetPatchMap),
        fSubstituteTargetVersion(substituteTargetVersion),
        fUseTargetFlags(useTargetFlags),
        fOnlyManualForTargetSubstituteOverlap(onlyManualForTargetSubstituteOverlap),
        fCurrentGeneNum(-1) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfWriter& mappedGxfFh,
                GxfWriter* unmappedGxfFh,
                ostream& mappingInfoFh,
                ostream* transcriptPslFh);
};

#endif
