/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "featureTree.hh"
#include "typeOps.hh"
#include <set>
class TransMap;
class PslMapping;
struct psl;
class PslCursor;
class AnnotationSet;
class BedMap;

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
    const BedMap* fTargetPatchMap; // location of patch regions in target genome
    const string fSubstituteTargetVersion;  // pass through targets when gene new gene doesn't map
    unsigned fUseTargetFlags;  // what targets to force.

    /* set of base ids (gene, transcript, havana) and gene names that have been
     * mapped.  Used to prevent output of target genes types that are not being
     * mapped (automatic genes), when type or source changes are already mapped.
     * N.B. Can't use AnnotationSet to track this, due to PAR mappings needing
     * to be mapped twice. */
    StringSet fMappedIdsNames;

    bool isSrcSeqInMapping(const GxfFeature* feature) const;
    bool isSrcSeqInMapping(const FeatureNode* featureNode) const;
    void recordMapped(const FeatureNode* featureNode);
    void recordGeneMapped(const FeatureNode* geneTree);
    bool checkMapped(const FeatureNode* featureNode) const;
    bool checkGeneMapped(const FeatureNode* geneTree) const ;

    ResultFeatureTrees processTranscript(const FeatureNode* transcriptTree,
                                         ostream* transcriptPslFh) const;
    ResultFeatureTreesVector processTranscripts(const FeatureNode* geneTree,
                                           ostream* transcriptPslFh) const;
    FeatureNode* findMatchingBoundingNode(const FeatureNodeVector& features,
                                          const FeatureNode* srcFeature) const;
    void copyMappingMetadata(const FeatureNode* origFeature,
                             FeatureNode* newFeature) const;
    void copyGeneMetadata(const FeatureNode* origGene,
                          FeatureNode* newGene) const;
    void forceToUnmapped(ResultFeatureTrees* mappedGene) const;
    void forceToUnmappedDueToRemapStatus(ResultFeatureTrees* mappedGene,
                                         RemapStatus remapStatus) const;
    void forceToUnmappedDueToTargetStatus(ResultFeatureTrees* mappedGene,
                                          TargetStatus targetStatus) const;
    bool hasMixedMappedSeqStrand(const ResultFeatureTrees* mappedGene) const;
    bool hasTargetStatusNonOverlap(const ResultFeatureTrees* mappedGene) const;
    bool hasExcessiveSizeChange(const ResultFeatureTrees* mappedGene) const;

    void setNumGeneMappings(FeatureNode* mappedGeneTree) const;
    bool shouldSubstituteTarget(const ResultFeatureTrees* mappedGene) const;
    FeatureNode* cloneTargetGene(const FeatureNode* srcGene) const;
    void substituteTarget(ResultFeatureTrees* mappedGene) const;
    void updateMappedGeneBounds(const FeatureNode* mappedTranscript,
                                string& seqid, string& strand,
                                int& start, int& end) const;
    FeatureNode* buildMappedGeneFeature(const FeatureNode* srcGeneTree,
                                        ResultFeatureTreesVector& mappedTranscripts) const;
    FeatureNode* buildUnmappedGeneFeature(const FeatureNode* srcGeneTree,
                                          ResultFeatureTreesVector& mappedTranscripts) const;
    ResultFeatureTrees buildGeneFeature(const FeatureNode* srcGeneTree,
                                        ResultFeatureTreesVector& mappedTranscripts) const;
    void saveMapped(ResultFeatureTrees& mappedGene,
                    AnnotationSet& mappedSet) const;
    void saveUnmapped(ResultFeatureTrees& mappedGene,
                      AnnotationSet& unmappedSet) const;
    void outputInfoHeader(ostream& mappingInfoFh) const;
    const GxfFeature* getTargetAnnotation(const FeatureNode* featureNode) const;
    const FeatureNode* getTargetAnnotationNode(const FeatureNode* featureNode) const;
    TargetStatus getTargetAnnotationStatus(const ResultFeatureTrees* mappedFeature) const;
    const string& getTargetAnnotationBiotype(const ResultFeatureTrees* mappedFeature) const;
    void outputFeatureInfo(const ResultFeatureTrees* mappedGene,
                           bool substituteTarget,
                           ostream& mappingInfoFh) const;
    void outputTranscriptInfo(const ResultFeatureTrees* mappedGene,
                              bool substituteTarget,
                              const FeatureNode* srcTranscript,
                              ostream& mappingInfoFh) const;
    void outputInfo(const ResultFeatureTrees* mappedGene,
                    ostream& mappingInfoFh) const;
    void processGeneLevelMapping(ResultFeatureTrees* mappedGene);
    void setGeneLevelMappingAttributes(ResultFeatureTrees* mappedGene);
    void mapGene(const FeatureNode* srcGeneTree,
                 AnnotationSet& mappedSet,
                 AnnotationSet& unmappedSet,
                 ostream& mappingInfoFh,
                 ostream* transcriptPslFh);
    RemapStatus getNoMapRemapStatus(const FeatureNode* geneTree) const;
    bool shouldMapGeneType(const FeatureNode* geneTree) const;
    bool inTargetPatchRegion(const FeatureNode* targetGene);
    bool checkTargetOverlappingMapped(const FeatureNode* targetGene,
                                      AnnotationSet& mappedSet);
    bool shouldIncludeTargetGene(const FeatureNode* geneTree,
                                 AnnotationSet& mappedSet);
    void copyTargetGene(const FeatureNode* targetGeneNode,
                        AnnotationSet& mappedSet,
                        ostream& mappingInfoFh);
    void copyTargetGenes(AnnotationSet& mappedSet,
                         ostream& mappingInfoFh);
    public:
    /* Constructor */
    GeneMapper(const AnnotationSet* srcAnnotations,
               const TransMap* genomeTransMap,
               const AnnotationSet* targetAnnotations,
               const BedMap* targetPatchMap,
               const string& substituteTargetVersion,
               unsigned useTargetFlags):
        fSrcAnnotations(srcAnnotations),
        fGenomeTransMap(genomeTransMap),
        fTargetAnnotations(targetAnnotations),
        fTargetPatchMap(targetPatchMap),
        fSubstituteTargetVersion(substituteTargetVersion),
        fUseTargetFlags(useTargetFlags) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfWriter& mappedGxfFh,
                GxfWriter& unmappedGxfFh,
                ostream& mappingInfoFh,
                ostream* transcriptPslFh);
};

#endif
