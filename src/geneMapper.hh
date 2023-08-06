/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "featureTree.hh"
#include "typeOps.hh"
#include <iostream>
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
    /* set of (baseId, chrom) or (name, chrom) that have been mapped. */
    class MappedIdSet: public set<pair<string, string> > {
        typedef pair<string, string> Key;
        public:
        void addBaseId(const string& fullid, const string& chrom) {
            insert(Key(getBaseId(fullid), chrom));
        }
        bool haveBaseId(const string& fullid, const string& chrom) const {
            return find(Key(getBaseId(fullid), chrom)) != end();
        }
        void removeBaseId(const string& fullid, const string& chrom) {
            erase(Key(getBaseId(fullid), chrom));
        }
        void addName(const string& name, const string& chrom) {
            insert(Key(name, chrom));
        }
        bool haveName(const string& name, const string& chrom) const {
            return find(Key(name, chrom)) != end();
        }
        void removeName(const string& name, const string& chrom) {
            erase(Key(name, chrom));
        }
        void dump(ostream& fh, const string& prefix="") const {
            for (auto key = cbegin(); key != cend(); key++) {
                fh << prefix << '(' << key->first << ',' << key->second << ')' << endl;
            }
        }
    };

    
    const AnnotationSet* fSrcAnnotations; // source annotations
    const TransMap* fGenomeTransMap;  // genomic mapping
    const AnnotationSet* fTargetAnnotations; // targeted genes/transcripts, maybe NULL
    const AnnotationSet* fPreviousMappedAnotations; // previous version
    const BedMap* fTargetPatchMap; // location of patch regions in target genome
    const string fSubstituteTargetVersion;  // pass through targets when gene new gene doesn't map
    unsigned fUseTargetFlags;  // what targets to force.
    bool fOnlyManualForTargetSubstituteOverlap;  // only check manual transcripts when checking target/map overlap

    /* set of base ids (gene, transcript, havana) and gene names that have been
     * mapped.  The key is "ident chrom" to handle PAR cases.
     * Used to prevent output of target genes types that are not being
     * mapped (automatic genes), when type or source changes are already mapped.
     */
    MappedIdSet fMappedIdsNames;
    
    int fCurrentGeneNum;  /* used by output info log to logically group features together,
                           * increments each time a gene is process */ 
    
    void outputInfoHeader(ostream& mappingInfoFh) const;
    void outputInfo(const string& recType,
                    const string& featType,
                    const FeatureNode* feature,
                    RemapStatus mappingStatus,
                    int mappingCount,
                    TargetStatus targetStatus,
                    ostream& mappingInfoFh) const;
    void outputSrcGeneInfo(const ResultFeatureTrees* mappedGene,
                           ostream& mappingInfoFh) const;
    void outputMappedGeneInfo(const ResultFeatureTrees* mappedGene,
                              ostream& mappingInfoFh) const;
    void outputUnmappedGeneInfo(const ResultFeatureTrees* mappedGene,
                                ostream& mappingInfoFh) const;
    void outputTargetGeneInfo(const ResultFeatureTrees* mappedGene,
                              const string& targetAction,
                              ostream& mappingInfoFh) const;
    string featureDesc(const FeatureNode* feature) const;
    bool isSrcSeqInMapping(const FeatureNode* feature) const;
    void debugRecordMapped(const FeatureNode* feature,
                           const string& desc,
                           const string& key = "") const;
    void recordGeneMapped(const FeatureNode* gene);
    void forgetGeneMapped(const FeatureNode* gene);
    void recordTranscriptMapped(const FeatureNode* transcript);
    void forgetTranscriptMapped(const FeatureNode* transcript);
    bool checkGeneMapped(const FeatureNode* gene) const ;
    bool checkTranscriptMapped(const FeatureNode* transcript) const;
    bool checkAllGeneTranscriptsMapped(const FeatureNode* gene) const;
    bool checkAnyGeneTranscriptsMapped(const FeatureNode* gene) const;
    ResultFeatureTrees processTranscript(const FeatureNode* transcript,
                                     ostream* transcriptPslFh);
    ResultFeatureTreesVector processTranscripts(const FeatureNode* gene,
                                            ostream* transcriptPslFh);
    FeatureNode* findMatchingBoundingFeature(const FeatureNodeVector& features,
                                             const FeatureNode* srcFeature) const;
    void copyMappingMetadata(const FeatureNode* origFeature,
                             FeatureNode* newFeature) const;
    void copyGeneMetadata(const FeatureNode* origGene,
                          FeatureNode* newGene) const;
    void forceToUnmapped(ResultFeatureTrees* mappedGene);
    void forceToUnmappedDueToRemapStatus(ResultFeatureTrees* mappedGene,
                                         RemapStatus remapStatus);
    void forceToUnmappedDueToTargetStatus(ResultFeatureTrees* mappedGene,
                                          TargetStatus targetStatus);
    bool hasMixedMappedSeqStrand(const ResultFeatureTrees* mappedGene) const;
    bool hasTargetStatusNonOverlap(const ResultFeatureTrees* mappedGene) const;
    bool hasExcessiveSizeChange(const ResultFeatureTrees* mappedGene) const;

    void setNumGeneMappings(FeatureNode* mappedGeneTree) const;
    bool checkForPathologicalGeneRename(const ResultFeatureTrees* mappedGene,
                                        const FeatureNode* targetGene) const;
    bool shouldSubstituteTarget(const ResultFeatureTrees* mappedGene) const;
    void substituteTarget(ResultFeatureTrees* mappedGene);
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
                    AnnotationSet& mappedSet);
    void saveUnmapped(ResultFeatureTrees& mappedGene,
                      AnnotationSet& unmappedSet);
    const FeatureNode* getTargetAnnotation(const FeatureNode* feature) const;
    TargetStatus getTargetAnnotationStatus(const ResultFeatureTrees* mappedFeature) const;
    const string& getTargetAnnotationBiotype(const ResultFeatureTrees* mappedFeature) const;
    void processGeneLevelMapping(ResultFeatureTrees* mappedGene);
    void setGeneLevelMappingAttributes(ResultFeatureTrees* mappedGene);
    void mapGene(const FeatureNode* srcGeneTree,
                 AnnotationSet& mappedSet,
                 AnnotationSet& unmappedSet,
                 FeatureTreePolish& featureTreePolish,
                 ostream& mappingInfoFh,
                 ostream* transcriptPslFh);
    void maybeMapGene(const FeatureNode* srcGeneTree,
                      AnnotationSet& mappedSet,
                      AnnotationSet& unmappedSet,
                      FeatureTreePolish& featureTreePolish,
                      ostream& mappingInfoFh,
                      ostream* transcriptPslFh);
    RemapStatus getNoMapRemapStatus(const FeatureNode* gene) const;
    bool shouldMapGeneType(const FeatureNode* gene) const;
    bool inTargetPatchRegion(const FeatureNode* targetGene);
    bool checkTargetOverlappingMapped(const FeatureNode* targetGene,
                                      AnnotationSet& mappedSet);
    bool checkForIncludeTargetSpecialCases(const FeatureNode* targetGene) const;
    bool shouldIncludeTargetGene(const FeatureNode* gene,
                                 AnnotationSet& mappedSet);
    void maybeCopyTargetGene(const FeatureNode* targetGene,
                             AnnotationSet& mappedSet,
                             ostream& mappingInfoFh);
    void copyTargetGene(const FeatureNode* targetGene,
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
                ostream& mappingInfoFh,
                ostream* transcriptPslFh);
};

#endif
