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
    SeqIdMap fMappedSeqRegionsWritten;  // mapped sequence ids that have been written

    /* check if a seqregion for seqid has been written, if so, return true,
     * otherwise record it and return false.  */
    bool checkRecordSeqRegionWritten(const string& seqid) {
        if (fMappedSeqRegionsWritten.find(seqid) == fMappedSeqRegionsWritten.end()) {
            fMappedSeqRegionsWritten.insert(pair<const string, bool>(seqid, true));
            return false;
        } else {
            return true;
        }
    }
    
    bool isSrcSeqInMapping(FeatureNode* featureNode) const;
    void processTranscript(FeatureNode* transcriptTree) const;
    void processTranscripts(FeatureNode* geneNode) const;
    bool haveMappedTranscripts(FeatureNode* geneNode) const;
    bool haveUnmappedTranscripts(FeatureNode* geneNode) const;
    void buildMappedGeneFeature(FeatureNode* geneNode,
                                bool srcSeqInMapping) const;
    void buildUnmappedGeneFeature(FeatureNode* geneNode,
                                  bool srcSeqInMapping) const;
    void outputMappedSeqRegionIfNeed(const GxfFeature* feature,
                                     ostream& mappedGxfFh);
    void outputMapped(FeatureNode* featureNode,
                      ostream& mappedGxfFh);
    void outputUnmapped(FeatureNode* featureNode,
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
