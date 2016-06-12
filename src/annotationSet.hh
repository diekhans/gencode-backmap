/*
 * map of location of target gene and transcripts. used in selecting
 * from multiple mappings
 */
#ifndef annotationSet_hh
#define annotationSet_hh

#include "gxf.hh"
#include <map>
#include <stdexcept>
#include "featureTree.hh"
struct genomeRangeTree;
class GenomeSizeMap;

/*
 * Locations in target genome of old transcripts, by base id
 */
class AnnotationSet {
    private:
    /* stored in range tree to link target location to
     * tree. */
    struct LocationLink {
        struct LocationLink *next;
        FeatureNode* featureNode;
    };
 
       
    // map of gene or transcripts to features. Keep up to two for PAR
    typedef map<const string, FeatureNodeVector> FeatureMap;
    typedef FeatureMap::iterator FeatureMapIter;
    typedef FeatureMap::const_iterator FeatureMapConstIter;

    // map by base id of genes and transcripts (not exons).
    FeatureMap fIdFeatureMap;

    // map by base id of exons
    FeatureMap fIdExonMap;

    // map by names of genes and transcripts
    FeatureMap fNameFeatureMap;

    // list of all gene features found
    FeatureNodeVector fGenes;

    // map of location to feature
    struct genomeRangeTree* fLocationMap;

    // mapped sequence ids that have been written
    StringSet fSeqRegionsWritten;

    // optional table of chromosome sequence sizes
    const GenomeSizeMap* fGenomeSizes;
    
    void addFeature(FeatureNode* featureNode);
    void processRecord(GxfParser *gxfParser,
                       GxfRecord* gxfRecord);
    void addLocationMap(FeatureNode* featureNode);
    void buildLocationMap();
    void freeLocationMap();
    bool isOverlappingGene(const FeatureNode* geneTree,
                           const FeatureNode* overlappingFeature,
                           float minSimilarity,
                           bool manualOnlyTranscripts);

    FeatureNode* getFeatureNodeByKey(const string& key,
                                     const FeatureMap& featureMap,
                                     const string& seqIdForParCheck) const;

    /* check if a seqregion for seqid has been written, if so, return true,
     * otherwise record it and return false.  */
    bool checkRecordSeqRegionWritten(const string& seqid) {
        if (fSeqRegionsWritten.find(seqid) == fSeqRegionsWritten.end()) {
            fSeqRegionsWritten.insert(seqid);
            return false;
        } else {
            return true;
        }
    }
    void outputSeqRegion(const string& seqId,
                         int size,
                         GxfWriter& gxfFh);
    void outputMappedSeqRegionIfNeed(const FeatureNode* geneTree,
                                     GxfWriter& mappedGxfFh);
    void outputFeature(const FeatureNode* featureNode,
                       GxfWriter& gxfFh) const;

    public:
    /* constructor, load gene and transcript objects from a GxF */
    AnnotationSet(const string& gxfFile,
                  const GenomeSizeMap* genomeSizes=NULL);

    /* constructor, empty set */
    AnnotationSet(const GenomeSizeMap* genomeSizes=NULL):
        fLocationMap(NULL),
        fGenomeSizes(genomeSizes) {
    }

    /* destructor */
    ~AnnotationSet();

    /* add a gene the maps */
    void addGene(FeatureNode* geneTree);

    /* get a gene or transcript node with same base id or NULL.  special
     * handling for PARs. */
    FeatureNode* getFeatureNodeById(const string& id,
                                    const string& seqIdForParCheck) const;

    /* get a gene or transcript node with same name or NULL.  special handling
     * for PARs. */
    FeatureNode* getFeatureNodeByName(const string& name,
                                      const string& seqIdForParCheck) const;

    /* get a target gene or transcript with same base or NULL.
     * special handling for PARs. */
    GxfFeature* getFeatureById(const string& id,
                               const string& seqIdForParCheck) const;

    /* get exon nodes by base id */
    FeatureNodeVector getExonNodesById(const string& exonId) const;

    /* find overlapping features */
    FeatureNodeVector findOverlappingFeatures(const string& seqid,
                                             int start,
                                             int end);
    
    /* find overlapping genes with minimum similarity at the transcript level */
    FeatureNodeVector findOverlappingGenes(const FeatureNode* geneTree,
                                           float minSimilarity,
                                           bool manualOnlyTranscripts);

    /* get list of all gene features */
    const FeatureNodeVector& getGenes() const {
        return fGenes;
    }

    /* print for debugging */
    void dump(ostream& fh) const;

    /* output genes */
    void write(GxfWriter& gxfFh);
};

#endif
