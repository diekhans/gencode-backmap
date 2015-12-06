/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "featureTree.hh"
#include <set>
class TransMap;
class PslMapping;
struct psl;
class PslCursor;
class TargetAnnotations;

/* class that maps a gene to the new assemble */
class GeneMapper {
    private:
    const TransMap* fGenomeTransMap;  // genomic mapping
    const TargetAnnotations* fTargetAnnotations; // targeted genes/transcripts, maybe NULL
    const string fSubstituteTargetVersion;  // pass through targets when gene new gene doesn't map
    bool fUseTargetForAutoGenes;  // don't map automatic genes
    bool fUseTargetForPseudoGenes;  // don't map pseudogenes

    typedef set<string> StringSet;
    StringSet fMappedSeqRegionsWritten;  // mapped sequence ids that have been written

    /* set of base ids (gene, transcript, havana) and gene names that have been
     * mapped.  Used to prevent output of target genes types that are not being
     * mapped (automatic genes), when type or source changes are already mapped */
    StringSet fMappedIdsNames;

    /* check if a seqregion for seqid has been written, if so, return true,
     * otherwise record it and return false.  */
    bool checkRecordSeqRegionWritten(const string& seqid) {
        if (fMappedSeqRegionsWritten.find(seqid) == fMappedSeqRegionsWritten.end()) {
            fMappedSeqRegionsWritten.insert(seqid);
            return false;
        } else {
            return true;
        }
    }

    bool isSrcSeqInMapping(const GxfFeature* feature) const;
    bool isSrcSeqInMapping(const FeatureNode* featureNode) const;
    void recordMapped(const FeatureNode* featureNode);
    void recordGeneMapped(const FeatureNode* geneTree);
    bool checkMapped(const FeatureNode* featureNode);
    bool checkGeneMapped(const FeatureNode* geneTree);
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
    void outputMappedSeqRegionIfNeed(const FeatureNode* geneTree,
                                     GxfWriter& mappedGxfFh);
    void outputFeature(const FeatureNode* featureNode,
                       GxfWriter& gxfFh) const;
    void outputSubstituted(const FeatureNode* featureNode,
                           GxfWriter& mappedGxfFh) const;
    void outputFeatures(const ResultFeatureTrees& mappedGene,
                        GxfWriter& mappedGxfFh,
                        GxfWriter& unmappedGxfFh);
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
    void mapGene(FeatureNode* srcGeneTree,
                 GxfFeature* geneFeature,
                 GxfWriter& mappedGxfFh,
                 GxfWriter& unmappedGxfFh,
                 ostream& mappingInfoFh,
                 ostream* transcriptPslFh);
    RemapStatus getNoMapRemapStatus(const FeatureNode* geneTree);
    bool shouldMapGeneType(const FeatureNode* geneTree);
    void copyTargetGene(const FeatureNode* targetGeneNode,
                        GxfWriter& mappedGxfFh,
                        ostream& mappingInfoFh);
    void copyTargetGenes(GxfWriter& mappedGxfFh,
                         ostream& mappingInfoFh);
    void processGene(GxfParser *gxfParser,
                     GxfFeature* geneFeature,
                     GxfWriter& mappedGxfFh,
                     GxfWriter& unmappedGxfFh,
                     ostream& mappingInfoFh,
                     ostream* transcriptPslFh);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord,
                       GxfWriter& mappedGxfFh,
                       GxfWriter& unmappedGxfFh,
                       ostream& mappingInfoFh,
                       ostream* transcriptPslFh);
    public:
    /* Constructor */
    GeneMapper(const TransMap* genomeTransMap,
               const TargetAnnotations* targetAnnotations,
               const string& substituteTargetVersion,
               bool useTargetForAutoGenes,
               bool useTargetForPseudoGenes):
        fGenomeTransMap(genomeTransMap),
        fTargetAnnotations(targetAnnotations),
        fSubstituteTargetVersion(substituteTargetVersion),
        fUseTargetForAutoGenes(useTargetForAutoGenes),
        fUseTargetForPseudoGenes(useTargetForPseudoGenes) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfParser *gxfParser,
                GxfWriter& mappedGxfFh,
                GxfWriter& unmappedGxfFh,
                ostream& mappingInfoFh,
                ostream* transcriptPslFh);
};

#endif
