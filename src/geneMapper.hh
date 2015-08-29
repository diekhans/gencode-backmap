/**
 * Handling mapping of a gene object
 */
#ifndef geneMapper_hh
#define geneMapper_hh
#include "gxf.hh"
#include "featureTree.hh"
#include <map>
class TransMap;
class PslMapping;
struct psl;
class PslCursor;

/* class that maps a gene to the new assemble */
class GeneMapper {
    private:
    const TransMap* fGenomeTransMap;  // genomic mapping
    
    typedef map<string, bool> SeqIdMap; // used for outputting GFF3 ##sequence-region metadata
    typedef pair<string, bool> SeqIdMapPair;
    SeqIdMap fMappedSeqRegionsWritten;  // mapped sequence ids that have been written

    /* check if a seqregion for seqid has been written, if so, return true,
     * otherwise record it and return false.  */
    bool checkRecordSeqRegionWritten(const string& seqid) {
        if (fMappedSeqRegionsWritten.find(seqid) == fMappedSeqRegionsWritten.end()) {
            fMappedSeqRegionsWritten.insert(SeqIdMapPair(seqid, true));
            return false;
        } else {
            return true;
        }
    }
    
    bool isSrcSeqInMapping(FeatureNode* featureNode) const;
    void processTranscript(FeatureNode* transcriptTree) const;
    void processTranscripts(FeatureNode* geneTree) const;
    void forceTranscriptToUnmapped(FeatureNode* featureNode,
                                   RemapStatus remapStatus) const;
    void forceTranscriptsToUnmapped(FeatureNode* geneTree,
                                    RemapStatus remapStatus) const;
    bool haveMappedTranscripts(FeatureNode* geneTree) const;
    bool haveUnmappedTranscripts(FeatureNode* geneTree) const;
    bool hasMixedMappedSeqStrand(FeatureNode* geneTree) const;
    int calcMappedGeneLength(FeatureNode* geneTree) const;
    bool hasExcessiveExpansion(FeatureNode* geneTree) const;

    void updateMappedGeneBounds(FeatureNode* transcriptTree,
                                string& seqid, string& strand,
                                int& start, int& end) const;
    void buildMappedGeneFeature(FeatureNode* geneTree,
                                bool srcSeqInMapping) const;
    void buildUnmappedGeneFeature(FeatureNode* geneTree,
                                  bool srcSeqInMapping) const;
    void buildGeneFeatures(FeatureNode* geneTree) const;
    void outputMappedSeqRegionIfNeed(const GxfFeature* feature,
                                     ostream& mappedGxfFh);
    void outputMapped(FeatureNode* featureNode,
                      ostream& mappedGxfFh);
    void outputUnmapped(FeatureNode* featureNode,
                        ostream& unmappedGxfFh);
    void output(FeatureNode* geneNode,
                ostream& mappedGxfFh,
                ostream& unmappedGxfFh);
    void processGene(GxfParser *gxfParser,
                     GxfFeature* geneFeature,
                     ostream& mappedGxfFh,
                     ostream& unmappedGxfFh);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord,
                       ostream& mappedGxfFh,
                       ostream& unmappedGxfFh);
    public:
    /* Constructor */
    GeneMapper(const TransMap* genomeTransMap):
        fGenomeTransMap(genomeTransMap) {
    }

    /* Map a GFF3/GTF */
    void mapGxf(GxfParser *gxfParser,
                ostream& mappedGxfFh,
                ostream& unmappedGxfFh);
};

#endif
